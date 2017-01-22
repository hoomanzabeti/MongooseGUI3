class LP:
    def findPosSupport():
        pass

    def __init__(self, opt, constraints, bounds):
        self.opt = opt
        self.constraints = constraints
        self.bounds = bounds
    def intro(Filename = 'trial.lp'):
        intro = 'Problem\n'
        intro += (Filename.replace('.lp','') + '\n')
        return intro
    def print_opt():
        if self.opt == Opt.MAXIMIZE:
            out = 'Maximize\n'
        elif self.opt == Opt.MINIMIZE:
            out = 'Minimize\n'
        elif self.opt == Opt.FEASIBLE:
            out = ''
        else:
            assert False
    def print_const():
        const = 'Subject\n'
        assert(False)
        #const += ''.join([name + ' ' + ])
    def print_bounds():
        bounds += 'Bounds\n'
        bounds += ''

    @staticmethod
    def findFreeLunch(N, Irrev, weight = [1], freeMetabs = []):
        m = len(N)
        n = len(N[0])

        sense = Sense.MAXIMIZE
        matrix = N
        space = 'row'
        if len(weight) == m:
            vectorObjective = [str(weight[i]) if weight[i] else 0 for i in range(m)]
        elif len(weight) == 1:
            vectorObjective = [1 for i in range(m)]
        else:
            raise Exception("Error: the weight vector is not of the right length!")
        



        #Opt
        if len(weight) == m:
            optimize_lin_comb = [('Y{1}'.format(str(i)), str(weight[i])) for i in range(m)]
        elif len(weight) == 1:
            optimize_lin_comb = [('Y{0}'.format(str(i)),1) for i in range(m)]
        else:
            raise Exception("Error: the weight vector is not of the right length!")
        opt = Opt(Opt.MAXIMIZE,optimize_lin_comb)
        

        #List of constraints
        constraints = []
        for i in range(m):
            lin_comb = [("X{0}".format(str(j)),str(N[i][j])) for j in range(n) if N[i][j]] + [(-1,'Y{0}'.format(str(i)))]
            bd_type = Constraint.EQ
            val = 0
            constraints += Constraint(lin_comb,bd_type,val)

        #List of bounds
        bounds = []
        lin_comb = [("X{0}".format(str(j)),1) for j in range(n) if j not in Irrev]
        bd_type = Constraint.NONE
        val = None
        bounds += Constraint(lin_comb,bd_type,val)
        assert(False)
        lin_comb = [()]

    @staticmethod
    def FBA_lp(N, growth, Exchange, allowed, limits = [1], forbidden = [], Negative = []):
        m = len(N)
        n = len(N[0])
        Forbidden = [x for x in (Exchange + forbidden) if x not in allowed]
        if not rec:
            Irrev = I
        else:
            Irrev = list(range(n))
        Rev = [x for x in range(n) if x not in Irrev]
        #Opt
        opt = Opt(Opt.MAXIMIZE,'V'+str(growth))
        #List of constraints
        constraints = []
        if Forbidden:
            lin_comb = [('V{0}'.format(str(i)),1) for i in Forbidden]
            bd_type = Constraint.EQ
            val = 0
            constraints.append(Constraint(lin_comb,bd_type,val))
        if len(limits) == 1:
            lin_comb = [('V{0}'.format(str(i)),1) for i in allowed]
            bd_type = Constraint.LEQ
            val = 1
            constraints.append(Constraint(lin_comb,bd_type,val))
        elif len(limits) == len(allowed):
            for i,xi in enumerate(allowed):
                lin_comb = [('V{0}'.format(xi),1)]
                bd_type = Constraint.LEQ
                val = xi
                constraints.append(Constraint(lin_comb,bd_type,val))
        else:
            raise Exception("Error: incompatible dimension of limits!")

        #List of bounds
        bounds = []
        if Negative:
            for i in Negative:
                lin_comb = [('V{0}'.format(str(i)),1)]
                bd_type = Constraint.LEQ
                val = 0
                bounds.append(Constraint(lin_comb,bd_type,val))
        for i in Rev:
            if i not in Negative and [_f for _f in [N[k][i] for k in range(m)] if _f]:
                lin_comb = [('V{0}'.format(str(i)),1)]
                bd_type = Constraint.NONE
                val = None
                bounds.append(Constraint(lin_comb,bd_type,val))
        self = LP.__init__(opt,constraints,bounds)

        return self




    def findPosSupport(N, support, weight = [1], Filename = 'trial.lp', Min = 0):
        pass

    def testCutSet(Cutset, N, Target, Filename = 'trial.lp', rec = True, I = []):
        m = len(N)
        n = len(N[0])
        intro = 'Problem\n'
        intro += (Filename.replace('.lp', '') + '\n')
        opt = 'Minimize\n\n'
        # note: we are only looking for a feasible vector, hence no objective function required!
        const = 'Subject\n'
        const += 'V' + str(Target) + ' = 1\n'
        const += ''.join([' + '.join([str(N[i][j]) + ' V' + str(j) for j in range(n) if N[i][j]]) + ' = 0\n' for i in range(m)])
        for react in Cutset:
            const += 'V' + str(react) + ' = 0\n'
        f = open(Filename, 'w')
        f.write(intro)
        f.write(opt)
        f.write(const)
        if not rec:
            bounds = 'Bounds\n'
            bounds += ''.join(['V' + str(j) + ' free\n' for j in range(n) if j not in I and [_f for _f in [N[i][j] for i in range(m)] if _f]])
            f.write(bounds)
        f.write('End\n')
        f.close()

class Opt():
    MAXIMIZE = 0
    MINIMIZE = 1
    FEASIBLE = 2
    def __init__(self,type,lin_comb = None):
        """lin_comb :: [(var_name,coeff)]
           type :: Opt.MAXIMIZE | Opt.MINIMIZE | Opt.FEASIBLE
        """
        self.type = type
        self.lin_comb = lin_comb

class Constraint():
    LEQ  = 0
    EQ   = 1
    GEQ  = 2
    NONE = 3
    def __init__(self, lin_comb, bd_type, val):
        """lin_comb :: [(var_name,coeff)]
           bd_type = Constraint.LEQ | Constraint.EQ | ...
           val = Int
        """
        self.lin_comb = lin_comb
        self.bd_type = bd_type
        self.val = val
    def printer(self):
        if self.type == Bound.LEQ:
            return ' <= {0}'.format(self.val)
        elif self.type == Bound.EQ:
            return ' = {0}'.format(self.val)
        elif self.type == Bound.GEQ:
            return ' >= {0}'.format(self.val)
        elif self.type == Bound.NONE:
            return ' free\n'
