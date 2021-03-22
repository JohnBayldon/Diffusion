
k_PA12 = 0.1
k_carbon = 1.0




#Rayleigh model for directioal fibers
def heat_capacity(k_matrix,k_fiber,V_f):
    k_comp = [0,0,0]
    k_comp[2] = k_matrix*(1+(k_fiber-k_matrix)/k_matrix)*V_f
    const_1 = (k_fiber+k_matrix)/(k_fiber-k_matrix)
    const_2 = 1/const_1
    k_comp[0] = k_matrix+1+2*V_f/(const_1-V_f + const_2*(0.30584*V_f**4+0.013363*V_f**8))
    return k_comp