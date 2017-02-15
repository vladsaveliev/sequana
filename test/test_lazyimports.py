# pytest Works locally but not at the package top level...


# This works but intervene with test_sphinxext somehow

def _test_lazy(mocker):
    import sys
    import sequana.lazy as lazy
    import sequana.lazyimports as lazyimports
    import imp
    import sphinx
    imp.reload(lazyimports)
    try:
        assert lazy.enabled() == False
    except Exception as err:
        raise(err)
    finally:
        del sys.modules['sphinx']
        imp.reload(lazyimports)
        assert lazy.enabled() == True
