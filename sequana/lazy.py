# Source inspiration and lazyimports.py taken from NIPY

from sequana.lazyimports import LazyImport

# lazy imports 
pylab = LazyImport('pylab')
numpy = LazyImport('numpy')
scipy = LazyImport('scipy')
scipy_stats = LazyImport('scipy.stats')
pandas = LazyImport('pandas')
reports = LazyImport('reports')
vcf = LazyImport("vcf")

def enabled():
    "Returns ``True`` if LazyImports are globally enabled"
    import nitime.lazyimports as l
    return not l.disable_lazy_imports
