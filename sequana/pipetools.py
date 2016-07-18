import os
import glob



class FixFailedProjects(object):
    """This is in progress but will search for
    directories that are not over yet and rebuild a multrun.sh file


    ::

        ff = FixFailedProjects()
        res = ff.create_new_starter()
        with open("multirun.sh") as fout:
            fout.write(res)

    """
    def __init__(self, where="."):
        
        self.pattern = "17 of 17"
        self.where = where
        self._template = """
        cd %(directory)s
        sh runme.sh &
        sleep 0.5
        echo Starting %(directory)s
        cd ..
        """

    def create_new_starter(self):
        count = 0
        txt = ""
        for directory in glob.glob("*"):
            if os.path.isdir(directory) is False:
                continue

            # else:
            filename = directory + "/run.err"; 
            if os.path.exists(filename) is False:
                txt += self._template % {"directory":directory}
                count += 1
            else:
                data = open(filename, "r").read()
                if self.pattern not in data:
                    txt += self._template % {"directory":directory}
                    count += 1
        return txt

