import matplotlib.pyplot as plt

import fem_1d

import os

# pure FEM:
#do_quad=False
#do_quad=True


import fem_1d

pp = fem_1d.parts(133)

mm = fem_1d.mesh(40)

pp.create_cos_f()
#pp.create_tophat_f()

#mm.onto_delta(pp.valueat_m)

mm.flip_volumes( pp )
mm.flip_assign( pp )

#plt.plot(pp.r , pp.f )

#mm.build_matrices()

#mm.onto_full( pp )

#plt.figure( figsize=[4,4] )
#plt.plot(mm.r , mm.f, 'o' )
#plt.show()

