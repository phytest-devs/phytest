import subprocess


def test_self_contained():
    testfile = 'examples/self_contained.py'
    p = subprocess.Popen(f"python {testfile}", stdout=subprocess.PIPE, shell=True)
    (output, err) = p.communicate()
    status = p.wait()
    assert status == 0
    assert '52 passed' in str(output)
