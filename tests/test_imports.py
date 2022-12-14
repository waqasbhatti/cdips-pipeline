def test_imports():

    modules = [
        'numpy',
        'astrobase',
        'astropy',
        'pandas',
        'astroquery',
        'imagesubphot',
        'aperturephot',
        'autoimagesub',
        'lcutils',
        'tessutils'
    ]

    for m in modules:

        dep_worked = True

        try:
            execstr =f"import {m} as foo"
            exec(execstr)
            dep_worked = True
        except Exception as e:
            print(e)
            dep_worked = False

        assert dep_worked, f"Failed at '{execstr}'"

    print("Success!")
    print("All import tests passed!")

if __name__ == "__main__":
    test_imports()
