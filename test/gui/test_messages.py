from sequana.gui import messages


def test_warning(qtbot):
    w = messages.WarningMessage("test")


def test_critical(qtbot):
    w = messages.CriticalMessage("test", details="test")




