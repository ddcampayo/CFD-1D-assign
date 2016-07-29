import matplotlib.pyplot as plt
import numpy as np

import fem_1d

#import os


pp = fem_1d.mesh(200)

mm = fem_1d.mesh(200)

#pp = fem_1d.parts(100)
#mm = fem_1d.parts(133)

mm.create_tophat_f()
pp.create_tophat_f()
#mm.create_cos_f()

#pp.create_cos_f()
#mm.onto_full( pp )

N = pp.r.size

v0 = 1.0
vel = v0 * np.ones( N )

dt = 0.1 / N

dr= dt * vel

T= 1 / v0

no_steps = int( T / dt) + 1

out = open('lagrangian.dat','w')

#pp.onto_full( mm )
#pp.onto_delta( mm )

pp.perturb(0.4)
#pp.randomize()

#mm.onto_delta( pp )

for step in xrange(0 , no_steps ) :
    time = step * dt

    #pp.onto_delta( mm )
    #pp.onto_lumped( mm )
    pp.onto_full( mm )
    #pp.onto_full_old( mm )

    pp.displace( dr )

    #mm.flip_volumes( pp )
    #mm.flip_assign( pp )

    #mm.onto_lumped( pp )
    #mm.onto_delta( pp )
    mm.onto_full( pp )
    #mm.onto_full_old( pp )
    
    out.write( "%f %f %f %f %f %f %f\n" %
               (time ,
                mm.int_f() , pp.int_f() , 
                mm.int_f2() , pp.int_f2() ,
                pp.fidelity(),mm.fidelity() ))

    if ( (step % (no_steps/10)) == 0 ):
        plt.plot( pp.r , pp.f )

plt.savefig('profiles.png')

out.close()


np.savetxt('part_final.dat' , zip( pp.r, pp.f  )  )
np.savetxt('mesh_final.dat' , zip( mm.r, mm.f  )  )

def move( r , dr):
    N = r.size
    r_new = np.zeros( N )

    for i in xrange(0,N):
        r_n = r[i] + dr[i]
        if (r_n > 1) :
            r_new[i] = r_n - 1
        elif (r_n < 0) :
            r_new[i] = r_n + 1
        else :
            r_new[i] = r_n

    return np.sort(r_new)
