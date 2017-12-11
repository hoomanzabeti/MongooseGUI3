def findPosSupport(N, support, weight = [1], Filename = 'trial.lp', Min = 0, restricted = True, option = 'row'):
    # This function finds the vector optimizing a given weight in the row/nullspace of N whose
    # support is restricted to a given set of entries; those entries must be non-negative!
    # Note: if the weight vector has a single component, it is automatically taken to be 1!
    # If a nonzero Min value is specified, the components in the support are at least Min.
    # If restricted = False, same but support is NOT restricted to the given set of entries.
    # print("NEW fps!")
    m, n = getSize(N)
    # intro = 'Problem\n'
    p = qsoptex.ExactProblem()
    # opt = 'Maximize\n'
    p.set_objective_sense(qsoptex.ObjectiveSense.MAXIMIZE)



    # if Min:
    #     if Min > 0:
    #         bounds += ''.join(['Y' + str(j) + ' >= ' + str(Min) + '\n' for j in support])
    #     else:
    #         bounds += ''.join([str(Min) + ' <= Y' + str(j) + ' <= ' + str(-Min) + '\n' for j in support])
    # else:
    #     bounds += ''.join(['Y' + str(j) + ' <= 1\n' for j in support])
    # if Cplex:
    #     const = const.replace('+ -', '-')

    if Min:
        if Min > 0:
            curLower, curUpper = Min, None
        else:
            curLower, curUpper = Min, -Min
    else:
        curLower, curUpper = 0, 1



    # Cplex is false
    # if Cplex:
    #     print()

    # if len(weight) == len(support):
    #     opt += ' + '.join([str(weight[i]) + ' Y' + str(support[i]) for i in range(len(support)) if weight[i]])
    # elif len(weight) == 1:   # equal weights, take all of them equal to 1
    #     opt += ' + '.join(['Y' + str(i) for i in support])
    # else:
    #     print('Error: the weight vector is not of the right length!')
    # opt += '\n'
    # const = 'Subject' + ' To' * int(Cplex) + '\n'
    variables = set([])

    # print(1 not in support)
    # exit()
    if len(weight) == len(support):
        for ind, item in enumerate(support):
            p.add_variable(name = 'Y' + str(item), objective = weight[ind], lower = curLower, upper = curUpper)
            # print(str(curLower) + ' <= ' + 'Y' + str(item) + ' <= ' + str(
            #     curUpper) + '    Objecitve coeficient\t' + str(weight[ind]))
            variables.add('Y' + str(item))
    elif len(weight) == 1:
        for ind, item in enumerate(support):
            p.add_variable(name = 'Y' + str(item), objective = 1, lower = curLower, upper = curUpper)
            # print(str(curUpper) + ' <= ' + 'Y' + str(item) + ' <= ' + str(curLower))
            variables.add('Y' + str(item))
    else:
        print('Error: the weight vector is not of the right length!')

    # if restricted:
    #     bounds += ''.join(['Y' + str(j) + ' = 0\n' for j in range(n) if j not in support])
    for ind in range(n):
        if ind not in support:
            p.add_variable(name = 'Y' + str(ind), objective = 0, lower = 0, upper = (0 if restricted else None))
            variables.add('Y' + str(ind))

    # if option == 'row':
    #     const += ''.join([' + '.join([str(N[i][j]) + ' X' + str(i) for i in range(m) if N[i][j]]) + ' + -1 Y' + str(j) + ' = 0\n' for j in range(n)])
    # else:
    #     const += ''.join([' + '.join([str(N[i][j]) + ' Y' + str(j) for j in range(n) if N[i][j]]) + ' = 0\n' for i in range(m)])
    # bounds = 'Bounds\n'
    # if option == 'row':
    #     bounds += ''.join(['X' + str(i) + ' free\n' for i in range(m) if [_f for _f in N[i] if _f]])

    if option == 'row':
        for ind in range(m):
            if [_f for _f in N[ind] if _f]:
                p.add_variable(name = 'X' + str(ind), objective = 0, lower = None, upper = None)
                variables.add('X' + str(ind))

        for j in range(n):
            curDict = {}
            for i in range(m):
                if N[i][j] != 0:
                    curDict.update({'X'+str(i) : N[i][j]})
                    variables.add('X'+str(i))
            curDict.update({'Y'+str(j): -1})
            variables.add('Y' + str(j))
            p.add_linear_constraint(qsoptex.ConstraintSense.EQUAL, curDict, rhs = 0)


    else: # option == 'null'
        for i in range(m):
            # print("NEW fps! (17)")
            curDict = {}
            for j in range(n):
                if N[i][j] != 0:
                    curDict.update({'Y'+str(j): N[i][j]})
                    curDict.update({'Y'+str(j): N[i][j]})
                    variables.add('Y' + str(j))
            p.add_linear_constraint(qsoptex.ConstraintSense.EQUAL, curDict, rhs=0)

    # opt += 'obj: '
    # f = open(Filename, 'w')
    # if not Cplex:
    #     f.write(intro)
    # f.write(opt)
    # f.write(const)
    # f.write(bounds)
    # f.write('End\n')
    # f.close()
    # return processFile(Filename, True)

    return processProblem(p, variables)
