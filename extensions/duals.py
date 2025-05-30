import os, sys, time
from pprint import pprint
from pyomo.environ import *
try:
    from pyomo.repn import generate_standard_repn
except ImportError:
    # this was called generate_canonical_repn before Pyomo 5.6
    from pyomo.repn import generate_canonical_repn as generate_standard_repn


def define_components(mod):

    # Make sure the model has dual and rc suffixes
    if not hasattr(mod, "dual"):
        mod.dual = Suffix(direction=Suffix.IMPORT)
    if not hasattr(mod, "rc"):
        mod.rc = Suffix(direction=Suffix.IMPORT)

def write_dual_costs(m):
    outputs_dir = m.options.outputs_dir

    # with open(os.path.join(outputs_dir, "producer_surplus{t}.csv".format(t=tag)), 'w') as f:
    #     for g, per in m.Max_Build_Potential:
    #         const = m.Max_Build_Potential[g, per]
    #         surplus = const.upper() * m.dual[const]
    #         if surplus != 0.0:
    #             f.write(','.join([const.name, str(surplus)]) + '\n')
    #     # import pdb; pdb.set_trace()
    #     for g, year in m.BuildGen:
    #         var = m.BuildGen[g, year]
    #         if var.ub is not None and var.ub > 0.0 and value(var) > 0.0 and var in m.rc and m.rc[var] != 0.0:
    #             surplus = var.ub * m.rc[var]
    #             f.write(','.join([var.name, str(surplus)]) + '\n')

    outfile = os.path.join(outputs_dir, "dual_costs.csv")
    dual_data = []
    start_time = time.time()
    print("Writing {} ... ".format(outfile), end=' ')

    def add_dual(m, const, lbound, ubound, duals, idx, prefix='', offset=0.0):
        if const in duals:
            dual = duals[const]
            if dual >= 0.0:
                direction = ">="
                bound = lbound
            else:
                direction = "<="
                bound = ubound

            # TODO: Maybe add something to undo the weighting on constraints defined for investment periods
            for i in idx:
                if i in m.TIMEPOINTS and i not in m.PERIODS:
                    dual = dual / m.bring_timepoint_costs_to_base_year[i]

            if bound is None:
                # Variable is unbounded; dual should be 0.0 or possibly a tiny non-zero value.
                if not (-1e-5 < dual < 1e-5):
                    raise ValueError("{} has no {} bound but has a non-zero dual value {}.".format(
                        const.name, "lower" if dual > 0 else "upper", dual))
            else:
                total_cost = dual * (bound + offset)
                if total_cost != 0.0:
                    dual_data.append((prefix+const.name, direction, (bound+offset), dual, total_cost))

    for comp in m.component_objects(ctype=Var):
        for idx in comp:
            var = comp[idx]
            if var.value is not None:  # ignore vars that weren't used in the model
                if var.is_integer() or var.is_binary():
                    # integrality constraint sets upper and lower bounds
                    add_dual(m, var, value(var), value(var), m.rc, idx, prefix='integer: ')
                else:
                    add_dual(m, var, var.lb, var.ub, m.rc, idx)
    for comp in m.component_objects(ctype=Constraint):
        for idx in comp:
            constr = comp[idx]
            if constr.active:
                offset = 0.0
                # cancel out any constants that were stored in the body instead of the bounds
                # (see https://groups.google.com/d/msg/pyomo-forum/-loinAh0Wx4/IIkxdfqxAQAJ)
                # (might be faster to do this once during model setup instead of every time)
                standard_constraint = generate_standard_repn(constr.body)
                if standard_constraint.constant is not None:
                    offset = -standard_constraint.constant                

                add_dual(m, constr, value(constr.lower), value(constr.upper), m.dual, idx, offset=offset)

    #dual_data.sort(key=lambda r: (not r[0].startswith('DR_Convex_'), r[3] >= 0)+r)

    with open(outfile, 'w') as f:
        f.write(','.join(['constraint', 'direction', 'bound', 'dual', 'total_cost']) + '\n')
        f.writelines(','.join(map(str, r)) + '\n' for r in dual_data)
        #f.writelines(','.join(map(str.replace(",", "."), r)) + '\n' for r in dual_data)
    print("time taken: {dur:.2f}s".format(dur=time.time()-start_time))

def post_solve(m, outdir):
    write_dual_costs(m)