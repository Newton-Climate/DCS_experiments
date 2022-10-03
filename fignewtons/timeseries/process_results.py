import h5py# declare physical constants


r = 8.314472; # universal gas constant joules / moles / K
N_a = 6.0221415e23; # molecules / moles
dz = 195017.0 # measurement pathlength in cm

def calc_vcd(p, T, dz, vmr_H2O=0.0):
    """
return vertical column density molec/cm^2
only dry vcd, no water 
"""
    rho_n = p*(1-vmr_H2O) / (r*T)*N_a/1.0e4
    vcd = dz*rho_n
    return vcd

def col2ppb(ch4, vcd, h2o):
    return 1e9 * ch4 / (vcd- h2o)

def col2ppm(co2, vcd, h2o):
    return 1e6 * co2 / (vcd- h2o)

def col2percent(h2o, vcd):
    return 1e2 * h2o / (vcd - h2o)



def read_data(datafile):
    with h5py.File(datafile, 'r') as f:
        read = lambda x: f[x][:]
        p, T = read('p'), read('T')
        h2o = read('h2o')
        vcd = calc_vcd(p, T, dz)
        co2, ch4 = read('co2'), read('ch4')

        out = {
            'CO2' : co2,
            'CO2_VMR' : col2ppm(co2, vcd, h2o),
            'CH4' : ch4,
            'CH4_VMR' : col2ppb(ch4, vcd, h2o),
            'H2O' : h2o,
            'H2O_VMR' : col2percent(h2o, vcd),
            'ch4_chi2' : read('ch4_chi2'),
            'co2_chi2' : read('co2_chi2'),
                        'p_ch4' : read('p_ch4'),
            'p' : p,
            'T' : T,
            'T_ch4' : read('T_ch4'),
                        'co2_measurement' : read('co2_measurement'),
            'ch4_measurement' : read('ch4_measurement'),
            'ch4_model' : read('ch4_model'),
            'co2_model' : read('co2_model')}
    return out

def read_OCO(datafile):
    with h5py.File(datafile, 'r') as f:
        read = lambda x: f[x][:]
        p, T = read('p'), read('T')
        h2o = read('h2o')
        vcd = calc_vcd(p, T, dz)
        co2 = read('co2')

        out = {
            'CO2' : co2,
            'CO2_VMR' : col2ppm(co2, vcd, h2o),
            'H2O' : h2o,
            'H2O_VMR' : col2percent(h2o, vcd),
            'co2_chi2' : read('co2_chi2'),
            'p' : read('p'),
            'T' : read('T'),
            'co2_measurement' : read('co2_measurement'),
            'co2_model' : read('co2_model')}
    return out

def read_obs(datafile):
    with h5py.File(datafile, 'r') as f:
        read = lambda x: f[x][:]

        out = {
            'CO2_VMR' : read('co2'),
            'H2O_VMR' : 1e2*read('h2o'),
            'CH4_VMR' : read('ch4'),
            't_pic' : read('t_pic'),
            'p' : read('p'),
            'T' : read('T'),
            't_pT' : read('t_pT')}
    return out
