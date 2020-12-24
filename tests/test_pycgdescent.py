import pytest
import pycgdescent as cg


def test_optimize_options():
    options = cg.OptimizeOptions()

    with pytest.raises(AttributeError):
        options.PrintLevel = 2

    options2 = options.replace(PrintLevel=2)
    assert options2.PrintLevel == 2

    print(options)
    print()
    print(options2)
    print()
    print(options2.pretty())


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        exec(sys.argv[1])
    else:
        from pytest import main
        main([__file__])

# vim: fdm=marker
