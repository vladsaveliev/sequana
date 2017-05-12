import re
import string
import subprocess
from collections import Counter, deque

from sequana.fasta import FastA
from sequana.lazy import pandas as pd
from sequana.lazy import numpy as np
from sequana.lazy import pylab


__all__ = ["DNA", "RNA", "Repeats", "Sequence"]


class Sequence(object):
    """Abstract base classe for other specialised sequences such as DNA.


    Sequenced is the base class for other classes such as :class:`DNA` and
    :class:`RNA`.

    ::

        from sequana import Sequence
        s = Sequence("ACGT")
        s.stats()
        s.get_complement()

    .. note:: You may use a Fasta file as input (see constructor)


    """
    def __init__(self, sequence, complement_in=b"ACGT", complement_out=b"TGCA",
                 letters="ACGT"):
        """.. rubric:: Constructor

        A sequence is just a string stored in the :attr:`sequence` attribute. It
        has properties related to the type of alphabet authorised.

        :param str sequence: May be a string of a Fasta File, in which case only
            the first sequence is used.
        :param complement_in:
        :param complement_out:
        :param letters: authorise letters. Used in :meth:`check` only.

        .. todo:: use counter only once as a property

        """
        if sequence.endswith(".fa") or sequence.endswith(".fasta"):
            fasta = FastA(sequence)
            sequence = fasta.next().sequence.upper()
        else: # assume correct string sequence
            pass

        self._data = sequence
        try:
            self._translate = string.maketrans(complement_in, complement_out)
        except:
            self._translate = bytes.maketrans(complement_in, complement_out)
        self._letters = letters

    def _get_sequence(self):
        return self._data
    sequence = property(_get_sequence)

    def get_complement(self):
        """Return complement """
        return self._data.translate(self._translate)

    def get_reverse_complement(self):
        """Return reverse complement """
        return self.get_complement()[::-1]

    def get_reverse(self):
        """Return reverse sequence"""
        return self._data[::-1]

    def complement(self):
        """Alias to :meth:`get_complement`"""
        self._data = self.get_complement()

    def reverse(self):
        """Alias to :meth:`get_reverse`"""
        self._data = self.get_reverse()

    def reverse_complement(self):
        """Alias to get_reverse_complement"""
        self._data = self.get_reverse_complement()

    def check(self):
        """Check that all letters are valid"""
        counter = Counter(self._data).keys()
        for key in counter:
            if key not in self._letters:
                raise ValueError("Found unexpected letter in the sequence (%s)" % key)

    def __len__(self):
        return len(self._data)

    def gc_content(self):
        """Return mean GC content"""
        c = Counter(self._data)
        ratio = (c['G'] + c['C']) / len(self.sequence)
        return ratio

    def stats(self):
        """Return basic stats about the number of letters"""
        from collections import Counter
        return Counter(self.sequence)

    def get_occurences(self, pattern, overlap=False):
        """Return position of the input pattern in the sequence

        ::

            >>> from sequana import Sequence
            >>> s = Sequence('ACGTTTTACGT')
            >>> s.get_occurences("ACGT")
            [0, 7]

        """
        if overlap is False:
            res = [m.start() for m in re.finditer(pattern, self.sequence)]
        elif overlap is True:
            res = [m.start() for m in re.finditer('(?=%s)'%pattern, self.sequence)]
        return res

        # reverse find-all without overlaps, you can combine positive and
        # negative lookahead into an expression like this:
        #res = [m.start() for m in re.finditer('(?=%s)(?!.{1,%d}%s)' % (search,
        #    len(pattern)-1, pattern), 'ttt')]


class DNA(Sequence):
    """Simple DNA class


        >>> d = DNA("ACGTTTT")
        >>> d.complement
        >>> d.reverse_complement

    """
    def __init__(self, sequence):
        super(DNA, self).__init__(sequence, complement_in=b"ACGT",
            complement_out=b"TGCA", letters="ACGTN")

        self._window       = None
        self._type_window  = None
        self._slide_window = None
        self._seq_right    = None
        self._dict_nuc     = dict_nuc = {'A' : 0, 'T' : 1, 'G' : 2, 'C' : 3}
        self._cumul        = None
        self._Xn           = None
        self._Yn           = None
        self._Zn           = None
        self._ignored_nuc  = None
        self._AT_skew_slide = None
        self._GC_skew_slide = None
        self._GC_content_slide = None
        self._AT_content_slide = None
        self._template_fft     = None
        self._c_fft            = None

    def _get_window(self):
        return self._window

    def _set_window(self,window):
        if (window > 0) & (window < 1):
            self._type_window = 'adapted to genome length : %.1f %% of total length' %(window*100)
            self._window = int(round(self.__len__() * window))

        elif (window >= 1) & ( window <= self.__len__()):
            self._type_window = 'fixed window length : %d' %(window)
            self._window = int(window)
        else:
            raise ValueError("Incorrect value for window: choose either float ]0,1]" +
                            " (fraction of genome) or integer [1,genome_length] (window size)")

        print("Computing GC skew")
        self._compute_skews()

        #### takes a too long time
        # print("Computing fft")
        # self._myfft_gc_skew(self._window)

    window = property(_get_window, _set_window)

    def _get_type_window(self):
        return self._type_window

    type_window = property(_get_type_window)

    def _init_sliding_window(self):
        """slide_window : deque of n first nucleotides (size of window)
        seq_right : deque of the rest of the genome, with the n-1 nucleotides at the end (circular DNA)
        """
        # split sequence_reference in two : sliding window and seq_right
        # for circular genome : add begining at the end of genome
        self._slide_window = deque(self.sequence[0:self._window])
        self._seq_right = deque( self.sequence[self._window:self.__len__()] + self.sequence[0:(self._window-1)] )

    def _init_list_results(self):
        # init IJ content and IJ skew
        IJ_content_res = np.empty((1,self.__len__()))
        IJ_content_res[:] = np.NAN
        IJ_skew_res = np.empty((1,self.__len__()))
        IJ_skew_res[:] = np.NAN
        return IJ_content_res, IJ_skew_res

    def _init_cumul_nuc(self):
        """Cumulative of nucleotide count along genome (init from first window)"""
        # ATGC (index stored in self._dict_nuc)
        cumul = np.zeros((4,(self.__len__()+self._window) ))
        for j in range(self._window):
            nuc = self.sequence[j]
            if nuc in self._dict_nuc:
                cumul[self._dict_nuc[nuc]][j] += 1
        self._cumul = cumul
    
    def _compute_skews(self):
        ### initialisation =  Calculating GC skew and AT skew for first window
        self._init_sliding_window()
        GC_content_slide, GC_skew_slide = self._init_list_results()
        AT_content_slide, AT_skew_slide = self._init_list_results()
        self._init_cumul_nuc()

        c = Counter(self._slide_window)
        dict_counts = {'G' : c['G'], 'C' : c['C'], 'A' : c['A'], 'T' : c['T']}
        i = 0

        # GC
        sumGC = float(dict_counts['G'] + dict_counts['C'])
        GC_content_slide[0][i] = sumGC
        if sumGC > 0:
            GC_skew_slide[0][i] = (dict_counts['G'] - dict_counts['C'])/sumGC
        # AT
        sumAT = float(dict_counts['A'] + dict_counts['T'])
        AT_content_slide[0][i] = sumAT
        if sumAT > 0:
            AT_skew_slide[0][i] = (dict_counts['A'] - dict_counts['T'])/sumAT

        ### Compute for all genome
        while(self._seq_right):
            out_nuc = self._slide_window.popleft()
            in_nuc = self._seq_right.popleft()
            self._slide_window.append(in_nuc)

            i += 1

            if i % 500000 == 0:
                print("%d / %d" % (i, self.__len__()))
            # if in and out are the same : do nothing, append same result
            if out_nuc != in_nuc:
                # remove out from counters
                if out_nuc in self._dict_nuc:
                    dict_counts[out_nuc] -= 1
                if in_nuc in self._dict_nuc:
                    dict_counts[in_nuc] += 1
                sumGC = float(dict_counts['G'] + dict_counts['C'])
                sumAT = float(dict_counts['A'] + dict_counts['T'])

            # fill results
            # GC
            GC_content_slide[0][i] = sumGC
            if sumGC > 0:
                GC_skew_slide[0][i] = (dict_counts['G'] - dict_counts['C'])/sumGC
            # AT
            AT_content_slide[0][i] = sumAT
            if sumAT > 0:
                AT_skew_slide[0][i] = (dict_counts['A'] - dict_counts['T'])/sumAT
            # cumul
            if in_nuc in self._dict_nuc:
                self._cumul[self._dict_nuc[in_nuc]][i+self._window-1] +=1

        self._GC_content_slide = GC_content_slide/float(self._window)
        self._AT_content_slide = AT_content_slide/float(self._window)
        self._cumul = np.delete(self._cumul, range(self.__len__(),self._cumul.shape[1]),1)
        self._cumul = np.cumsum(self._cumul,axis=1)

        ### save result for Z curve
        self._Xn = list((self._cumul[self._dict_nuc['A']] + self._cumul[self._dict_nuc['G']]) -\
         (self._cumul[self._dict_nuc['C']] + self._cumul[self._dict_nuc['T']]))

        self._Yn = list((self._cumul[self._dict_nuc['A']] + self._cumul[self._dict_nuc['C']]) -\
         (self._cumul[self._dict_nuc['G']] + self._cumul[self._dict_nuc['T']]))

        self._Zn = list((self._cumul[self._dict_nuc['A']] + self._cumul[self._dict_nuc['T']]) -\
         (self._cumul[self._dict_nuc['C']] + self._cumul[self._dict_nuc['G']]))

        self._AT_skew_slide = AT_skew_slide
        self._GC_skew_slide = GC_skew_slide

        ### check proportion of ignored nucleotides
        GC_content_total = (self._cumul[self._dict_nuc['G']][-1] + self._cumul[self._dict_nuc['C']][-1]) / float(self.__len__())
        AT_content_total = (self._cumul[self._dict_nuc['A']][-1] + self._cumul[self._dict_nuc['T']][-1]) / float(self.__len__())
        self._ignored_nuc = 1.0 - GC_content_total - AT_content_total

    def _get_AT_skew(self):
        if self._AT_skew_slide is None:
            raise AttributeError("Please set a valid window to compute skew")
        else:
            return self._AT_skew_slide

    AT_skew = property(_get_AT_skew)

    def _get_GC_skew(self):
        if self._GC_skew_slide is None:
            raise AttributeError("Please set a valid window to compute skew")
        else:
            return self._GC_skew_slide

    GC_skew = property(_get_GC_skew)


    def _create_template_fft(self,M=1000):
        M_3 =  int(M/3)
        W = [-0.5] * M_3 + list(np.linspace(-0.5,0.5,M-2*M_3)) + [0.5] * M_3
        return list(W * np.hanning(M))

    def _myfft_gc_skew(self,M):
        """
        x : GC_skew vector (list)
        param N: length of the GC skew vector
        param M: length of the template
        param A: amplitude between positive and negative GC skew vector

        """
        x = self._GC_skew_slide[0]

        N = len(x)
        template =  self._create_template_fft(M) + [0] * (N-M)
        template/=pylab.norm(template)
        print("1")
        c = np.fft.fft(x)
        print("1bis")

        c = abs(np.fft.ifft(
                        np.fft.fft(x) * pylab.conj(np.fft.fft(template))
                        )**2)/pylab.norm(x)/pylab.norm(template)
        print("2")
        # shift the SNR vector by the template length so that the peak is at the END of the template
        c = np.roll(c, M//2)
        print("3")
        self._template_fft = template
        self._c_fft        = c*2./N

    def plot_all_skews(self,figsize=(10, 12), fontsize=16, alpha=0.5):
        if self._window is None:
            raise AttributeError("Please set a valid window to compute skew")

        # create figure
        # fig, axarr = pylab.subplots(10,1, sharex=True, figsize=figsize)
        fig, axarr = pylab.subplots(9,1, sharex=True, figsize=figsize)

        main_title = "Window size = %d (%.0f %% of genome )\n\
        GC content = %.0f %%, AT content = %.0f %%, ignored = %.0f %%" \
        % (self._window, self._window*100/self.__len__(), 
            self.gc_content()*100, (1-self.gc_content())*100, self._ignored_nuc*100)

        pylab.suptitle(main_title, fontsize=fontsize)

        # GC skew
        axarr[0].set_title("GC skew (blue) - Cumulative sum (red)")
        axarr[0].plot(list(self._GC_skew_slide[0]),'b-',alpha=alpha)
        axarr[0].set_ylabel("(G -C) / (G + C)")

        axarr[1].plot(list(np.cumsum(self._GC_skew_slide[0])),'r-',alpha=alpha)
        axarr[1].set_ylabel("(G -C) / (G + C)")
        
        # AT skew
        axarr[2].set_title("AT skew (blue) - Cumulative sum (red)")
        axarr[2].plot(list(self._AT_skew_slide[0]),'b-',alpha=alpha)
        axarr[2].set_ylabel("(A -T) / (A + T)")

        axarr[3].plot(list(np.cumsum(self._AT_skew_slide[0])),'r-',alpha=alpha)
        axarr[3].set_ylabel("(A -T) / (A + T)")

        # Xn
        axarr[4].set_title("Cumulative RY skew (Purine - Pyrimidine)")
        axarr[4].plot(self._Xn,'g-',alpha=alpha)
        axarr[4].set_ylabel("(A + G) - (C + T)")

        # Yn
        axarr[5].set_title("Cumulative MK skew (Amino - Keto)")
        axarr[5].plot(self._Yn,'g-',alpha=alpha)
        axarr[5].set_ylabel("(A + C) - (G + T)")

        # Zn
        axarr[6].set_title("Cumulative H-bond skew (Weak H-bond - Strong H-bond)")
        axarr[6].plot(self._Zn,'g-',alpha=alpha)
        axarr[6].set_ylabel("(A + T) - (G + C)")

        # GC content
        axarr[7].set_title("GC content")
        axarr[7].plot(list(self._GC_content_slide[0]),'k-',alpha=alpha)
        axarr[7].set_ylabel("GC")

        # AT content
        axarr[8].set_title("AT content")
        axarr[8].plot(list(self._AT_content_slide[0]),'k-',alpha=alpha)
        axarr[8].set_ylabel("AT")

        # # FFT
        # axarr[9].set_title("FFT")
        # axarr[9].plot(list(self._c_fft),'g-',alpha=alpha)
        # axarr[9].set_ylabel("FFT")

        fig.tight_layout()
        fig.subplots_adjust(top=0.88)





class RNA(Sequence):
    """Simple RNA class


        >>> d = RNA("ACGUUUU")
        >>> d.complement
        >>> d.reverse_complement

    """
    def __init__(self, sequence):
        super(RNA, self).__init__(sequence, complement_in=b"ACGU",
            complement_out=b"UGCA", letters="ACGUN")


class Repeats(object):
    """Class for finding repeats in DNA or RNA linear sequences.

    Computation is performed each time the :attr:`threshold` is set 
    to a new value.

    .. plot::
        :include-source:

        from sequana import sequana_data, Repeats
        rr = Repeats(sequana_data("measles.fa"))
        rr.threshold = 4
        rr.hist_length_repeats()

    .. note:: Works with shustring package from Bioconda (April 2017)
    .. todo:: use a specific sequence (right now it is the first one)

    """
    def __init__(self, filename_fasta, merge=False):
        """.. rubric:: Constructor

        Input must be a fasta file with valid DNA or RNA characters

        :param str filename_fasta: a Fasta file, only the first
            sequence is used !
        :param int threshold: Minimal length of repeat to output
        """
        self._threshold                       = None
        self._df_shustring                    = None
        self._header                          = None
        self._length                          = None
        self._longest_shustring               = None
        self._begin_end_repeat_position       = None
        self._begin_end_repeat_position_merge = None
        self._filename_fasta                  = filename_fasta
        self._previous_thr                    = None
        self._list_len_repeats                = None
        if not isinstance(merge, bool):
            raise TypeError("do_merge must be boolean")
        self._do_merge                        = merge

    def _get_header(self):
        """get first line of fasta (needed in input shustring)
        and replace spaces by underscores
        """
        if self._header is None:
            # -n is to specify the DNA nature of the sequence
            first_line = subprocess.check_output(["head","-n","1",self._filename_fasta])
            first_line = first_line.decode('utf8')
            first_line = first_line.replace("\n","").replace(" ","_")
            self._header = first_line
        return self._header
    header = property(_get_header)

    def _get_shustrings_length(self):
        """Return dataframe with shortest unique substring length at each position
        shortest unique substrings are unique in the sequence and its complement
        Uses shustring tool"""
        if self._df_shustring is None:
            # read fasta
            task_read           = subprocess.Popen(["cat", self._filename_fasta],
                stdout=subprocess.PIPE)
            # replace spaces of fasta
            task_replace_spaces = subprocess.Popen(["sed","s/ /_/g"],
                stdin=task_read.stdout, stdout=subprocess.PIPE)
            # shustring command
            task_shus           = subprocess.Popen(['shustring','-r','-q','-l',self.header],
                stdin=task_replace_spaces.stdout, stdout=subprocess.PIPE)

            # read stdout line by line and append to list
            list_df = []
            for line in task_shus.stdout:
                list_df.append(line.decode('utf8').replace("\n",'').split("\t"))
                #df=pd.concat([df, line])
            task_shus.wait()

            # convert to dataframe
            df = pd.DataFrame(list_df[2:len(list_df)])
            self._df_shustring = df.astype(int)
            self._df_shustring.columns = ["position","shustring_length"]

            # get input sequence length and longest shustring in the first line
            self._length = int(list_df[0][1])
            self._longest_shustring = int(list_df[0][3].split("<=")[2])

        return self._df_shustring
    df_shustring = property(_get_shustrings_length)

    def _get_genome_length(self):
        if self._df_shustring is None:
            self._get_shustrings_length()
        return self._length
    length = property(_get_genome_length)

    def _get_longest_shustring(self):
        if self._df_shustring is None:
            self._get_shustrings_length()
        return self._longest_shustring
    longest_shustring = property(_get_longest_shustring)

    def _find_begin_end_repeats(self,force=False):
        """Returns position of repeats longer than threshold
        as an ordered list
        """
        if self.df_shustring is None:
            self._get_shustrings_length()

        if self._threshold is None:
            #print("No threshold : please set minimul length of repeats to output")
            raise ValueError("threshold : please set threshold (minimum length of repeats to output)")

        # if there is no result yet, or the threshold has changed
        if (self._begin_end_repeat_position is None) | (self.threshold != self._previous_thr) | force:
            nb_row = self.df_shustring.shape[0]
            i = 0
            step_repeat_seq = []
            be_repeats = []
            e = 0

            # use list because faster
            list_len_shus = list(self.df_shustring.loc[:,"shustring_length"])

            while(i < nb_row):
                # begining of repeat
                if (list_len_shus[i] > self.threshold):
                    b = i
                    # compute new end of repeat
                    len_repeat = list_len_shus[i]
                    e = b + len_repeat
                    # save (b,e)
                    be_repeats.append((b,e))
                    # update i
                    i = e-1
                i +=1

            self._begin_end_repeat_position = be_repeats

        self._get_merge_repeats()

    def _get_be_repeats(self):
        self._find_begin_end_repeats()
        return self._begin_end_repeat_position

    begin_end_repeat_position = property(_get_be_repeats)

    def _set_threshold(self,value):
        if value < 0:
            raise ValueError("Threshold must be positive integer")
        if value != self._threshold:
            self._previous_thr = self._threshold
        self._threshold = value
        self._find_begin_end_repeats()
        self._list_len_repeats = [tup[1]-tup[0] for tup in self._begin_end_repeat_position]

    def _get_threshold(self):
        return self._threshold

    threshold = property(_get_threshold, _set_threshold)

    def _get_list_len_repeats(self):
        if self._list_len_repeats is None:
            raise UserWarning("Please set threshold (minimum length of repeats to output)")
        return self._list_len_repeats

    list_len_repeats = property(_get_list_len_repeats)

    def _get_merge_repeats(self):
        if self._do_merge:
            # if there are repeats, merge repeats that are fragmented
            if len(self._begin_end_repeat_position) > 0:
                prev_tup = self._begin_end_repeat_position[0]
                b = prev_tup[0]
                begin_end_repeat_position_merge = []
                for i in range(1,len(self._begin_end_repeat_position)):
                    tup = self._begin_end_repeat_position[i]

                    if tup[0] == prev_tup[1]:
                        # concat
                        e = tup[1]
                        prev_tup = tup
                        if i == (len(self._begin_end_repeat_position) -1):
                            # last tup : append to result
                            begin_end_repeat_position_merge.append((b,e))

                    else:
                        # real end of repeat : append result and update b, e
                        e = prev_tup[1]
                        begin_end_repeat_position_merge.append((b,e))
                        prev_tup = tup
                        b = prev_tup[0]
                        if i == (len(self._begin_end_repeat_position) -1):
                            # last tup : append to result
                            begin_end_repeat_position_merge.append(tup)

                self._begin_end_repeat_position = begin_end_repeat_position_merge


    def _get_do_merge(self):
        return self._do_merge

    def _set_do_merge(self, do_merge):
        if not isinstance(do_merge, bool):
            raise TypeError("do_merge must be boolean")
        # if different
        if do_merge != self._do_merge:
            self._do_merge = do_merge
            if self._do_merge:
                # did not merge before, merge now
                if self._begin_end_repeat_position is None:
                    self._find_begin_end_repeats()
            else:
                # data is already merged : need to compute again to un-merge
                self._find_begin_end_repeats(force=True)
    do_merge = property(_get_do_merge,_set_do_merge)

    def hist_length_repeats(self, bins=None, alpha=0.5, hold=False,
            fontsize=12, grid=True, label="Repeat length",
            xlabel="Repeat length", ylabel="#"):
        """Plots histogram of the repeat lengths


        """
        # check that user has set a threshold
        if self._list_len_repeats is None:
            self._get_list_len_repeats()

        if bins is None:
            bins = range(max(0, self.threshold -1), max(self._list_len_repeats)+2)

        if hold is False:
            pylab.clf()
        pylab.hist(self._list_len_repeats, alpha=alpha, label=label, bins=bins)
        pylab.legend()
        pylab.xlabel(xlabel, fontsize=fontsize)
        pylab.ylabel(ylabel, fontsize=fontsize)
        if grid is True:
            pylab.grid(True)

