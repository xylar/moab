import sys
import traceback

class colors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def test_driver(test_list):
    ret_val = 0
    for test in test_list:
        try:
            test()
        except:
            print(colors.FAIL + "FAIL" + colors.ENDC + ": " + test.__name__)
            ret_val += 1
            traceback.print_exc()
        else:
            print(colors.OKGREEN + "PASS" + colors.ENDC + ": " + test.__name__)
    sys.exit(ret_val)

def CHECK_EQ(actual_value, expected_value):
    err_msg = "Expected value: {} Actual value: {}"
    err_msg = err_msg.format(expected_value, actual_value)
    assert(expected_value == actual_value, err_msg)

def CHECK_NOT_EQ(actual_value, expected_value):
    err_msg = "Expected value: not {} Actual value: {}"
    err_msg = err_msg.format(expected_value, actual_value)
    assert(expected_value != actual_value, err_msg)
