import sys

def arg_passed(arg_name):
    # Takes abbreviations into account
    possible_matches = []
    for passed_arg in sys.argv[1:]:
        if(passed_arg.strip()==arg_name[:len(passed_arg.strip())]):
            possible_matches.append(passed_arg)
    argname_passed = len(possible_matches) == 1
    return argname_passed
