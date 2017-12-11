import sys
import shutil
import unittest
import qsoptex
from ModelProcessing import processFile

class Tester(unittest.TestCase):
    def test_all(self):
        # print('HERE!')
        f1 = open('tests/test01_tester_01.txt', 'r')
        f2 = open('tests/test01_judge_01.txt', 'r')
        content_tester = ''
        for line in f1.readlines():
            content_tester += line
        content_judge = ''
        for line in f2.readlines():
            content_judge += line

        self.assertMultiLineEqual(content_judge, content_tester)
        f1.close()
        f2.close()

        # # print('HERE!')
        # f1 = open('tests/test01_tester_02.txt', 'r')
        # f2 = open('tests/test01_judge_02.txt', 'r')
        # content_tester = ''
        # for line in f1.readlines():
        #     content_tester + line
        # content_judge = ''
        # for line in f2.readlines():
        #     content_judge + line
        #
        # self.assertMultiLineEqual(content_judge, content_tester)
        # f1.close()
        # f2.close()


if '-r' in sys.argv:
    p = qsoptex.ExactProblem()
    p.set_objective_sense(qsoptex.ObjectiveSense.MAXIMIZE)
    p.add_variable(name="x", objective=1, lower=0, upper=7)
    p.add_variable(name="y", objective=1, lower=0, upper=9)
    p.add_linear_constraint(qsoptex.ConstraintSense.LESS,
                            {'x': 1, 'y': 1}, rhs=15)
    status = p.solve()
    if status == qsoptex.SolutionStatus.OPTIMAL:
        value = p.get_objective_value()
        print("Qsoptex Python API's answer: ")
        print(value)
        print("x: " + str(p.get_value('x')))
        print("y: " + str(p.get_value('y')))


    intro = 'Problem\ntest\n'
    opt = 'Maximize\nobj: 1 x + 1 y\n'
    const = 'Subject\n1 x + 1 y <= 15\n'
    bounds = 'Bounds\nx <= 7\ny <= 9\n'
    f = open('test.lp', 'w')
    f.write(intro)
    f.write(opt)
    f.write(const)
    f.write(bounds)
    f.write('End\n')
    f.close()
    print("Executive qsoptex answer: ")
    processFile('test.lp', True)
    exit()

elif '-j' in sys.argv:
    open('mode.txt', 'w').write('0')
elif '-c' in sys.argv:
    unittest.main()
    Tester.test_all()
    exit()

elif '-t' in sys.argv:
    ind = sys.argv.index('-t')
    if len(sys.argv) != ind+1 and sys.argv[ind+1].isdigit():
        open('mode.txt', 'w',).write(sys.argv[ind+1])
    else:
        open('mode.txt', 'w', ).write('-1')

import tests.Test01

if '-j' in sys.argv:
    shutil.move('Reduction.txt', 'tests/test01_judge_01.txt')
    shutil.move('ReductionFull.txt', 'tests/test01_judge_02.txt')
elif '-t' in sys.argv:
    shutil.move('Reduction.txt', 'tests/test01_tester_01.txt')
    shutil.move('ReductionFull.txt', 'tests/test01_tester_02.txt')

