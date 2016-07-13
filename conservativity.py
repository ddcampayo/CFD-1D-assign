import matplotlib.pyplot as plt

import fem_1d

#import os

# pure FEM:
#do_quad=False
#do_quad=True

pp = fem_1d.mesh(100)

mm = fem_1d.mesh(133)

pp.create_cos_f()
pp.check_cos_f()
pp.int_f()
pp.int_f2()

#mm.flip_volumes( pp )
#mm.flip_assign( pp )


mm.onto_full( pp )
mm.check_cos_f()
mm.int_f()
mm.int_f2()

plt.plot( pp.r , pp.f )
plot(mm.r, mm.f, 'o')

pp.onto_full( mm )
pp.check_cos_f()
pp.int_f()
pp.int_f2()
plot( pp.r , pp.f , '*' )
