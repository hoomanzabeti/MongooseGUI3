import sys
import json
from fractions import *

variables = set([])


def calc(formula, solution):
    global variables
    parts = formula.split()
    sum = 0
    i = 0
    sign = 1
    while i < len(parts):
        coefficient = Fraction(parts[i])
        variable = parts[i+1]
        if i+2 < len(parts) and parts[i+2] == '-':
            sign = -1
        else:
            sign = 1
        sum += coefficient * sign * solution[variable]
        variables.add(variable)
        i += 3
    return sum


def reversed_signed(sign):
    if sign == '=':
        return '='
    elif sign == '<=':
        return '>='
    elif sign == '<':
        return '>'
    elif sign == '>=':
        return '<='
    elif sign == '>':
        return '<'
    elif sign == 'free':
        return 'free'
    else:
        print('Error sgn: ' + sign)
        exit()


def __main__():
    problem_file = sys.argv[1]
    solution_file = sys.argv[2]

    p = open(problem_file, 'r')
    s = open(solution_file, 'r')
    solution_value = Fraction(s.readline())
    solution = {}
    solution_string = list(s.readline())
    frac = False
    for i in range(len(solution_string)):
        if solution_string[i] == '(':
            frac = True
        if frac and solution_string[i] == ',':
            solution_string[i] = '/'
            frac = False
    solution_string = "".join(solution_string)
    solution_string = solution_string.replace('Fraction(', '')
    solution_string = solution_string.replace(')', '')
    solution_string = solution_string.replace('/ ', '/')
    for entity in solution_string.split(","):
        parts = entity.split(":")
        variable = parts[0].split("'")[1]
        value = parts[1]
        solution.update({variable: Fraction(value)})
    s.close()
    problem = ''.join([line for line in p.readlines()])
    p.close()

    objective = problem.partition('obj: ')[2].partition('Subject')[0]
    consts = problem.partition('Subject')[2].partition('Bounds')[0]
    bounds = problem.partition('Bounds')[2].partition('End')[0]
    copy_of_bounds = bounds
    for bound in bounds.split('\n'):
        parts = bound.split()
        if len(parts) == 3 and (parts[1] == '<=' or parts[1] == '<') and Fraction(parts[2]) > 0:
            copy_of_bounds = copy_of_bounds.replace(bound, "0 <= " + bound)
        if len(parts) == 3 and (parts[1] == '>=' or parts[1] == '>') and Fraction(parts[2]) > 0:
            copy_of_bounds = copy_of_bounds.replace(bound, "0 >= " + bound)
    bounds = copy_of_bounds
    for bound in bounds.split('\n'):
        parts = bound.split()
        if len(parts) > 3:
            bound1 = parts[2] + " " + reversed_signed(parts[1]) + " " + parts[0]
            bound2 = parts[2] + " " + parts[3] + " " + parts[4]
            copy_of_bounds = copy_of_bounds.replace(bound, bound1+"\n"+bound2)
    bounds = copy_of_bounds
    # print(bounds)
    # exit()

    if calc(objective, solution) == solution_value:
        print('Correct answer for objective value')
    else:
        print('Value of objective is ' + str(calc(objective, solution[1])) + ', not ' + str(solution[0]))
        exit()
    for const in consts.split('\n'):

        if '=' in const:
            parts = const.partition("=")
            if calc(parts[0], solution) == Fraction(parts[2]):
                print('satisfied constraint')
            else:
                print('unsatisfied constraint:')
                print(const)
                print(calc(parts[0], solution))
                exit()

        if '<' in const:
            parts = const.partition("<")
            if calc(parts[0], solution[1]) < Fraction(parts[2]):
                print('satisfied constraint')
            else:
                print('unsatisfied constraint:')
                print(const)
                print(calc(parts[0], solution))
                exit()

        if '>' in const:
            parts = const.partition(">")
            if calc(parts[0], solution) > Fraction(parts[2]):
                print('satisfied constraint')
            else:
                print('unsatisfied constraint:')
                print(const)
                print(calc(parts[0], solution))
                exit()

        if '<=' in const:
            parts = const.partition("<=")
            if calc(parts[0], solution) <= Fraction(parts[2]):
                print('satisfied constraint')
            else:
                print('unsatisfied constraint:')
                print(const)
                print(calc(parts[0], solution))
                exit()

        if '>=' in const:
            parts = const.partition(">=")
            if calc(parts[0], solution) >= Fraction(parts[2]):
                print('satisfied constraint')
            else:
                print('unsatisfied constraint:')
                print(const)
                print(calc(parts[0], solution))
                exit()

    for bound in bounds.split('\n'):
        if bound == "":
            continue
        parts = bound.split()
        variable = parts[0]
        sign = parts[1]
        if sign == '=':
            if solution[variable] == Fraction(parts[2]):
                print('Bound satisfied')
            else:
                print('Bound did not satisfy: ')
                print(solution[variable])
                print(bound)
                exit()
        elif sign == '<':
            if solution[variable] < Fraction(parts[2]):
                print('Bound satisfied')
            else:
                print('Bound did not satisfy: ')
                print(solution[variable])
                print(bound)
                exit()
        elif sign == '>':
            if solution[variable] > Fraction(parts[2]):
                print('Bound satisfied')
            else:
                print('Bound did not satisfy: ')
                print(solution[variable])
                print(bound)
                exit()
        elif sign == '<=':
            if solution[variable] <= Fraction(parts[2]):
                print('Bound satisfied')
            else:
                print('Bound did not satisfy: ')
                print(solution[variable])
                print(bound)
                exit()
        elif sign == '>=':
            if solution[variable] >= Fraction(parts[2]):
                print('Bound satisfied')
            else:
                print('Bound did not satisfy: ')
                print(solution[variable])
                print(bound)
                exit()
        elif sign == 'free':
            print('Bound satisfied')
        else:
            print('Bound error: ' + bound)
            exit()
        if variable in variables:
            variables.remove(variable)

    for variable in variables:
        if solution[variable] < 0:
            print('Default bound violated:')
            print(variable + ": " + str(solution[variable]))
            exit()
        else:
            print('Default bound satisfied')

__main__()