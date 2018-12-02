#! /usr/bin/env/python3

def save_object(gnc_struct, name):
    with open(name, 'w') as fout:
        l = ['alpha','V_fps','beta','Nz','alpha_dot','beta_dot','del_el',
             'del_ail','del_rud','t1/2_lon','t1/2_lat']
        s = ','.join(l) + '\n'
        fout.write(s)
        for i in range(gnc_struct.len_count):
            s = ''
            s += "{:1.6f}".format(gnc_struct.alpha_rad[i]) + ','
            s += "{:6.0f}".format(gnc_struct.V_fps[i]) + ','
            s += "{:1.6f}".format(gnc_struct.beta_rad[i]) + ','

