def Reading():                  #Function for reading raw data and input parameters form CSV files
    with open('CSV/Raw_Data.csv','r') as Raw_Data_File:
        raw_data=csv.reader(Raw_Data_File)
        raw_data.__next__()
        for line in raw_data:
            Time.append(float(line[0]))
            Pressure.append(float(line[1]))
            Rate.append(float(line[2]))
    with open('CSV/Parameters.csv', 'r') as Parameters_File:
        para = csv.reader(Parameters_File)
        para.__next__()
        for line in para:
            Input.append(float(line[1]))

def Extract():   #Extraxt data
    global P_Inj
    global EndInj
    global Inj_Total

    P_Inj=Pressure[0]           #Pressure at first of Injection
    i=0
    x=0
    while i < len(Time)-1:      #Cumulative calculation
        i+=1
        x=Rate[i]*(Time[i]-Time[i-1])*60+x
        Cum.append(x)
        if Rate[i]==0:
            EndInj=Time[i-1]      #Shut-in time
            Inj_Total=Cum[-1]
            break
    dt_Shut.append(Time[i]-EndInj)
    Prs_Shut.append(Pressure[i])
    Last_Prs=Prs_Shut[0]
    while i < len(Time)-1:      #Shut-in data
        i+=1
        if Last_Prs-Pressure[i]>Res_Pres_Inc:
            dt_Shut.append(Time[i]-EndInj)
            Prs_Shut.append(Pressure[i])
            Last_Prs=Pressure[i]

def Compliance():
    global Slope
    global P_res
    global min_dpdg
    global T_max_dpdg
    global G_max_dpdg
    global P_max_dpdg
    global T_min_dpdg
    global G_min_dpdg
    global P_min_dpdg
    global ISIP
    global dpdt_1_2
    global dpdt_1

    i=0
    while i<=len(dt_Shut)-1:
        dp.append(Prs_Shut[0]-Prs_Shut[i])
        dtd.append(dt_Shut[i]/EndInj)
        G_time.append(4 / 3.14159 * 4 / 3 * ((1 + dtd[i]) ** 1.5 - dtd[i] ** 1.5 - 1))
        log_dt.append(math.log10(dt_Shut[i]))
        t_1_2.append(dt_Shut[i] ** (-1 / 2))
        t_1.append(dt_Shut[i] ** (-1))
        i += 1
    i=1
    dp_dg_1.append([])
    gdp_dg_1.append([])
    dp_dt.append([])
    tdp_dt.append([])
    while i<=len(dt_Shut)-2:
        dp_dg_1.append(-(Prs_Shut[i + 1] - Prs_Shut[i - 1]) / (G_time[i + 1] - G_time[i - 1]))
        gdp_dg_1.append(G_time[i] * dp_dg_1[i])
        dp_dt.append(-(Prs_Shut[i + 1] - Prs_Shut[i - 1]) / (dt_Shut[i + 1] - dt_Shut[i - 1]))
        tdp_dt.append(dt_Shut[i] * dp_dt[i])
        log_tdpdt.append(math.log10(tdp_dt[i]))
        i+=1
    Draft=5
    Slope = (log_tdpdt[-Draft]-log_tdpdt[-Draft*2])/(log_dt[-Draft]-log_dt[-Draft*2])
    dpdt_1_2 = (Prs_Shut[-Draft] - Prs_Shut[-1]) / (t_1_2[-Draft] - t_1_2[-1])
    dpdt_1 = (Prs_Shut[-Draft] - Prs_Shut[-1]) / (t_1[-Draft] - t_1[-1])
    if Slope>-0.75:
        P_res = -dpdt_1_2 * t_1_2[-Draft] + Prs_Shut[-Draft]
    else:
        P_res = -dpdt_1 * t_1[-Draft] + Prs_Shut[-Draft]
#''''''''''''''''Finding ISP
    j=1
    while dp_dg_1[-j] < dp_dg_1[-j - 1]:
        j += 1
    i = len(dt_Shut) - j-1
    P_max_dpdg = Prs_Shut[i]
    G_max_dpdg = G_time[i]
    T_max_dpdg = dt_Shut[i]
    while dp_dg_1[-j]>dp_dg_1[-j-1]:
        j+=1
    i=len(dt_Shut)-j-1
    min_dpdg=dp_dg_1[i]
    P_min_dpdg = Prs_Shut[i]
    G_min_dpdg = G_time[i]
    T_min_dpdg = dt_Shut[i]
    ISIP = P_min_dpdg + min_dpdg * G_min_dpdg
    i=1
    P_ISIP.append([])
    Flow.append([])
    while i<len(dt_Shut)-1:
        P_ISIP.append(Prs_Shut[i]-ISIP)
        Flow.append(WSC*dp_dt[i]/60)
        i+=1

def Tangent():
    global ShminT

    i=1
    while gdp_dg_1[-i-1]>gdp_dg_1[-i]:
        i+=1
    Flag=True
    while Flag:
        if gdp_dg_1[-i]>G_time[i]*(gdp_dg_1[-i]-gdp_dg_1[-i-1])/(G_time[-i]-G_time[-i-1]):
            Flag=False
        i+=1
    ShminT=Prs_Shut [-i-2]
    return -i-1

def h_function():
    global h_Shut
    global h_Peak

    i=0
    while i<len(dt_Shut):
        if dt_Shut[i]<T_min_dpdg:
            Eff_Prs.append(P_min_dpdg + (G_min_dpdg - G_time[i]) * min_dpdg)
        else:
            Eff_Prs.append(Prs_Shut[i])
        Sum=float(0)
        j=1
        while j < len(dt_Shut):
            if dt_Shut[j] < dt_Shut[i]:
                x=(Eff_Prs[j] - Eff_Prs[j - 1]) * (dt_Shut[i] - dt_Shut[j]) ** 0.5
                Sum += x
                j += 1
            else:
                break
        Matrix.append(Sum)
        i += 1
    i=0
    while i<len(dt_Shut):
        Term1.append((Eff_Prs[0] - P_res) * (EndInj * 0.5 + dt_Shut[i]) ** 0.5)
        h_fun.append(Term1[i] + Matrix[i])
        if dt_Shut[i]==T_max_dpdg:
            h_Peak=h_fun[i]
        i+=1
    h_Shut=h_fun[0]
    i=1
    Stiff.append([])
    while i<len(dt_Shut)-1:
        Stiff.append(-(Eff_Prs[i+1] - Eff_Prs[i-1]) / (h_fun[i+1] - h_fun[i-1]))
        i+=1

def Results():
    global strain
    global sf_PKN
    global sf_RDL
    global CL_PKN
    global CL_RDL
    global length_PKN
    global length_RDL
    global K_PKN
    global K_RDL
    global Mass_PKN
    global Mass_RDL
    global sf_PKN_T
    global sf_RDL_T
    global CL_PKN_T
    global CL_RDL_T
    global length_PKN_T
    global length_RDL_T
    global K_PKN_T
    global K_RDL_T
    global Mass_PKN_T
    global Mass_RDL_T
    global perm_RDL
    global area_Perm
    global A_K_pst
    global rad
    global radial_pst
    global perm_RDL_pst
    global area
    global area_pst
    global P_K_R
    global length
    global length_pst
    global perm_PKN
    global perm_PKN_pst
    global ShminC
    global G_K_C
    global G_K_T
    global H_K
    global P_K

    l_min=1.0
    l_max=1000.0
    inc=0.01

    strain = Young / (1 - Poisson ** 2)
    ShminC = ISIP-100

    # PKN Compliance
    sf_PKN = 2 * strain / (3.14159 * Frac_h)
    length_PKN =np.arange(l_min,l_max,inc) #452.122965031195
    CL_PKN = min_dpdg * (2 * WSC / 3.14159 * 5.615 + length_PKN * Frac_h ** 2 / strain) / (
            length_PKN * Frac_h * (EndInj * 60) ** 0.5)
    K_PKN = 1E+15 * (
            (CL_PKN / 3.2808 / 60 ** 0.5) / (ISIP - P_res) * 145.04) ** 2 * 3.14159 * Fluid_Vis * 0.000000001 / (
                    Poro * (Fluid_Com + Comp) * 145.04)
    Mass_PKN = abs(Inj_Total * 5.615 - ((ISIP - ShminC) * WSC + 4 * CL_PKN * length_PKN * Frac_h * (
            EndInj * 60 / 2) ** 0.5 + length_PKN * Frac_h * (ISIP - ShminC) / sf_PKN))
    i=0
    while Mass_PKN[i]>Mass_PKN[i+1]: #abs(int(Mass_PKN[i]*Resol)/Resol-0.5)>0:
        i+=1
    length_PKN[0]=length_PKN[i]
    CL_PKN[0]= CL_PKN[i]
    K_PKN[0]=K_PKN[i]
    Mass_PKN[0]=Mass_PKN[i]

    # Radial Compliance
    length_RDL = np.arange(l_min,l_max,inc) #61.0922808
    sf_RDL=3*3.14159*strain/16/length_RDL
    CL_RDL=min_dpdg*(WSC*5.615+length_RDL**3*16/3/strain)/(2*3.14159*length_RDL**2*(EndInj*60)**0.5)
    K_RDL=1E+15*((CL_RDL/3.2808/60**0.5)/(ISIP-P_res)* 145.04)**2*3.14159*Fluid_Vis*0.000000001/(Poro*(Comp+Fluid_Com)*145.04)
    Mass_RDL=abs(Inj_Total*5.615-((ISIP-ShminC)*WSC+4*CL_RDL*3.14159*length_RDL**2*(EndInj*60/2)**0.5+length_RDL**2*3.14159*(ISIP-ShminC)/sf_RDL))
    i=0
    while(Mass_RDL[i]>Mass_RDL[i+1]):
        i+=1
    sf_RDL[0]=sf_RDL[i]
    length_RDL[0] = length_RDL[i]
    CL_RDL[0] = CL_RDL[i]
    K_RDL[0] = K_RDL[i]
    Mass_RDL[0] = Mass_RDL[i]

    # PKN Tangent
    sf_PKN_T = 2 * strain / (3.14159 * Frac_h)
    length_PKN_T =np.arange(l_min,l_max,inc) #452.122965031195
    CL_PKN_T = min_dpdg * (2 * WSC / 3.14159 * 5.615 + length_PKN_T * Frac_h ** 2 / strain) / (
            length_PKN_T * Frac_h * (EndInj * 60) ** 0.5)
    K_PKN_T = 1E+15 * (
            (CL_PKN_T / 3.2808 / 60 ** 0.5) / (ISIP - P_res) * 145.04) ** 2 * 3.14159 * Fluid_Vis * 0.000000001 / (
                    Poro * (Fluid_Com + Comp) * 145.04)
    Mass_PKN_T = abs(Inj_Total * 5.615 - ((ISIP - ShminT) * WSC + 4 * CL_PKN_T * length_PKN_T * Frac_h * (
            EndInj * 60 / 2) ** 0.5 + length_PKN_T * Frac_h * (ISIP - ShminT) / sf_PKN_T))
    i=0
    while Mass_PKN_T[i]>Mass_PKN_T[i+1]: #abs(int(Mass_PKN[i]*Resol)/Resol-0.5)>0:
        i+=1
    length_PKN_T[0]=length_PKN_T[i]
    CL_PKN_T[0]= CL_PKN_T[i]
    K_PKN_T[0]=K_PKN_T[i]
    Mass_PKN_T[0]=Mass_PKN_T[i]

    # Radial Tangent
    length_RDL_T = np.arange(l_min,l_max,inc) #61.0922808
    sf_RDL_T=3*3.14159*strain/16/length_RDL_T
    CL_RDL_T=min_dpdg*(WSC*5.615+length_RDL_T**3*16/3/strain)/(2*3.14159*length_RDL_T**2*(EndInj*60)**0.5)
    K_RDL_T=1E+15*((CL_RDL_T/3.2808/60**0.5)/(ISIP-P_res)* 145.04)**2*3.14159*Fluid_Vis*0.000000001/(Poro*(Comp+Fluid_Com)*145.04)
    Mass_RDL_T=abs(Inj_Total*5.615-((ISIP-ShminT)*WSC+4*CL_RDL_T*3.14159*length_RDL_T**2*(EndInj*60/2)**0.5+length_RDL_T**2*3.14159*(ISIP-ShminT)/sf_RDL_T))
    i=0
    while(Mass_RDL_T[i]>Mass_RDL_T[i+1]):
        i+=1
    sf_RDL_T[0] = sf_RDL_T[i]
    length_RDL_T[0] = length_RDL_T[i]
    CL_RDL_T[0] = CL_RDL_T[i]
    K_RDL_T[0] = K_RDL_T[i]
    Mass_RDL_T[0] = Mass_RDL_T[i]

    # h-function
    # radial
    area_Perm = 3.2808 ** 3 * (0.9 * (Inj_Total / 6.29 - WSC / 6.29 * (P_max_dpdg - P_Inj))) / (
            4 * h_Peak / 145.04 * 60 * (
            Poro * (Comp + Fluid_Com) * 145.04 / (3.14159 * Fluid_Vis * 0.000000001)) ** 0.5)
    rad = ((3 * 3.14159 * strain / 145.04 * (
            Inj_Total / 6.29 - WSC / 6.29 * 145.04 * (ISIP - P_Inj) / 145.04 - 4 * area_Perm / 3.2808 ** 3 * (
            ISIP - P_res) / 145.04 * (EndInj / 2 * 3600) ** 0.5 * (
                    Poro * (Comp + Fluid_Com) * 145.04 / (3.14159 * Fluid_Vis * 0.000000001)) ** 0.5)) / (
                   16 * 3.14159 * (ISIP - ShminC) / 145.04)) ** (1 / 3) * 3.2808
    perm_RDL = (area_Perm / (3.14159 * rad ** 2)) ** 2 / 3.2808 ** 2 * 1E+15

    # PKN
    area_Meter = area_Perm / 3.2808 ** 3
    area = 3.2808 ** 2 * (Inj_Total / 6.29 - (ISIP - P_Inj) * WSC / 6.29 - area_Meter * 4 * h_Shut / 145.04 * 60 * (
            Poro * (Comp + Fluid_Com) * 145.04 / 3.14159 / (
            Fluid_Vis * 0.000000001)) ** 0.5) * sf_PKN/ 145.04 * 3.2808 / ((ISIP - ShminC) / 145.04)
    length = area / Frac_h
    perm_PKN = (area_Perm / 3.2808 ** 3 / (length / 3.2808 * Frac_h / 3.2808)) ** 2 * 1E+15

    # Postclosure
    A_K_pst = 3.2808 ** 3 * (dpdt_1_2 / 145.04 * 60) ** -1 * (
            Inj_Total / 6.29 - WSC / 6.29 * 145.04 * (P_res - P_Inj) / 145.04) / 2 * (
                      Fluid_Vis * 0.000000001 / (3.14159 * (Comp + Fluid_Com) * 145.04 * Poro)) ** 0.5
    # radial
    radial_pst = ((3 * 3.14159 * strain / 145.04 * (
            Inj_Total / 6.29 - WSC / 6.29 * 145.04 * (ISIP - P_Inj) / 145.04 - 4 * area_Perm / 3.2808 ** 3 * (
            ISIP - P_res) / 145.04 * (EndInj / 2 * 3600) ** 0.5 * (
                    Poro * (Comp + Fluid_Com) * 145.04 / (3.14159 * Fluid_Vis * 0.000000001)) ** 0.5)) / (
                          16 * 3.14159 * (ISIP - ShminC) / 145.04)) ** (1 / 3) * 3.2808
    perm_RDL_pst = (A_K_pst / (3.14159 * radial_pst ** 2)) ** 2 / 3.2808 ** 2 * 1E+15
    # PKN
    area_pst = 3.2808 ** 2 * (
            Inj_Total / 6.29 - (ISIP - P_Inj) * WSC / 6.29 - (A_K_pst / 3.2808 ** 3) * 4 * h_Shut / 145.04 * 60 * (
            Poro * (Comp + Fluid_Com) * 145.04 / 3.14159 / (
            Fluid_Vis * 0.000000001)) ** 0.5) * sf_PKN / 145.04 * 3.2808 / ((ISIP - ShminC) / 145.04)
    length_pst = area_pst / Frac_h
    perm_PKN_pst = (A_K_pst / 3.2808 ** 3 / (length_pst / 3.2808 * Frac_h / 3.2808)) ** 2 * 1E+15

    #Final value of permeability
    if Frac_h<0:
        G_K_C=K_RDL[0]
        G_K_T=K_RDL_T[0]
        H_K=perm_RDL
        P_K=perm_RDL_pst
    else:
        G_K_C = K_PKN[0]
        G_K_T = K_PKN_T[0]
        H_K = perm_PKN
        P_K = perm_PKN_pst

def Eport_CSV():
    #Export Estmations values
    with open('CSV/Estimation.csv','w', newline="") as Estimation_csv:
        EST_CSV=csv.writer(Estimation_csv)
        EST_CSV.writerows([['Item', 'Values'],
                           ['Injection Duration (hr)', EndInj],
                           ['Pressure at start of Injection', P_Inj],
                           ['Minimum dP/dG (psi)', min_dpdg],
                           ['G-time at minimum dP/dG (unitless)', G_min_dpdg],
                           ['Pressure at minimum dP/dG (psi)', P_min_dpdg],
                           ['Effective ISIP (psi)', ISIP],
                           ['G-function value at peak dP/dG', G_max_dpdg],
                           ['Pressure at peak dP/dG', P_max_dpdg],
                           ['h-function at shut-in (psi-hrs^(1/2))', h_Shut],
                           ['h-function at peak dP/dG (psi-hrs^(1/2))', h_Peak],
                           ['Reservoir Pressure', P_res],
                           ['Minimum principal stress (psi) (Compliance)', ShminC],
                           ['Minimum principal stress (psi) (Tangent)   ', ShminT],
                           ['dP/d(t^(-1/2)) (if present) (psi/hours^(-1/2))', dpdt_1_2],
                           ['dP/d(t^-1) (if present)', dpdt_1],
                           ['Impulse', Slope]
                           ])

    #Export PKN Geometry results from Compliance method
    with open('CSV/PKN_Compliance.csv','w',newline="") as Results_CSV:
        Results=csv.writer(Results_CSV)
        Results.writerows([['Item', 'Values'],
                           ['Sf (psi/ft)',sf_PKN],
                           ['Length (ft)', length_PKN[0]],
                           ['CL (ft/min^.5)', "{:.2e}".format(CL_PKN[0])],
                           ['Permeability (md)', "{:.2e}".format(K_PKN[0])],
                           ['Mass balance residual', Mass_PKN[0]]
                           ])


    # Export RDL Geometry results from Compliance method
    with open('CSV/RDL_Compliance.csv', 'w', newline="") as Results_CSV:
        Results = csv.writer(Results_CSV)
        Results.writerows([['Item', 'Values'],
                           ['Sf (psi/ft)', sf_RDL[0]],
                           ['Radius (ft)', length_RDL[0]],
                           ['CL (ft/min^.5)', "{:.2e}".format(CL_RDL[0])],
                           ['Permeability (md)', "{:.2e}".format(K_RDL[0])],
                           ['Mass balance residual', Mass_RDL[0]]
                           ])

    # Export PKN Geometry results from Tangent method
    with open('CSV/PKN_Tangent.csv', 'w', newline="") as Results_CSV:
        Results = csv.writer(Results_CSV)
        Results.writerows([['Item', 'Values'],
                           ['Sf (psi/ft)', sf_PKN_T],
                           ['Length (ft)', length_PKN_T[0]],
                           ['CL (ft/min^.5)', "{:.2e}".format(CL_PKN_T[0])],
                           ['Permeability (md)', "{:.2e}".format(K_PKN_T[0])],
                           ['Mass balance residual', Mass_PKN_T[0]]
                           ])

    # Export RDL Geometry results from Tangent method
    with open('CSV/RDL_Tangent.csv', 'w', newline="") as Results_CSV:
        Results = csv.writer(Results_CSV)
        Results.writerows([['Item', 'Values'],
                           ['Sf (psi/ft)', sf_RDL_T[0]],
                           ['Radius (ft)', length_RDL_T[0]],
                           ['CL (ft/min^.5)', "{:.2e}".format(CL_RDL_T[0])],
                           ['Permeability (md)', "{:.2e}".format(K_RDL_T[0])],
                           ['Mass balance residual', Mass_RDL_T[0]]
                           ])

    # Export RDL Geometry results from h-function method
    with open('CSV/RDL_h_function.csv', 'w', newline="") as Results_CSV:
        Results = csv.writer(Results_CSV)
        Results.writerows([['Item', 'Values'],
                           ['Area*k^(1/2) (ft^3)', area_Perm],
                           ['Radius (ft)', rad],
                           ['Permeability (md)', "{:.2e}".format(perm_RDL)]
                           ])

    # Export PKN Geometry results from h-function method
    with open('CSV/PKN_h_function.csv', 'w', newline="") as Results_CSV:
        Results = csv.writer(Results_CSV)
        Results.writerows([['Area*k^(1/2) (ft^3)',area_Perm],
                           ['Area estimate (ft^2)', area],
                           ['Length (ft)', length],
                           ['Permeability (md)', "{:.2e}".format(perm_PKN)]
                           ])

    # Export RDL Geometry results from postclosure method
    with open('CSV/RDL_Postclosure.csv', 'w', newline="") as Results_CSV:
        Results = csv.writer(Results_CSV)
        Results.writerows([['Area*k^(1/2) (ft^3)',A_K_pst],
                           ['Radius (ft)',radial_pst],
                           ['Permeability (md)',"{:.2e}".format(perm_RDL_pst)]
                           ])

    # Export PKN Geometry results from from postclosure method
    with open('CSV/PKN_Postclosure.csv', 'w', newline="") as Results_CSV:
        Results = csv.writer(Results_CSV)
        Results.writerows([['Area*k^(1/2) (ft^3)',A_K_pst],
                           ['Area estimate (ft^2)',area_pst],
                           ['Length (ft)',length_pst],
                           ['Permeability (md)',"{:.2e}".format(perm_PKN_pst)]
                           ])

    # Export calculated permeability via all methods
    with open('CSV/Permeability_Values.csv', 'w', newline="") as Results_CSV:
        Results = csv.writer(Results_CSV)
        Results.writerows([['G-function permeability (md) (Compliance)',"{:.2e}".format(G_K_C) ],
                           ['G-function permeability (md) (Tangent)   ',"{:.2e}".format(G_K_T)],
                           ['h-function permeability (md)',"{:.2e}".format(H_K) ],
                           ['Postclosure linear permeability (md)',"{:.2e}".format(P_K) ]
                           ])

    # Export Vectors
    with open('CSV/Vectors.csv', 'w', newline="") as Results_CSV:
        Results = csv.writer(Results_CSV)
        Results.writerow(['Time', 'Pressure', 'G_Time', 'dp/dG-1','G*dpdG-1'])
        i=0
        while i<len(dp_dg_1):
            Exp_Lst = []
            Exp_Lst.append(dt_Shut[i])
            Exp_Lst.append(Prs_Shut[i])
            Exp_Lst.append(G_time[i])
            Exp_Lst.append(dp_dg_1[i])
            Exp_Lst.append(gdp_dg_1[i])
            i+=1
            Results.writerow(Exp_Lst)


#Defining arra
Time=[]
Pressure=[]
Rate=[]
Cum=[]
dt_Shut=[]
Prs_Shut=[]
dp=[]
dtd=[]
G_time=[]
dp_dg_1=[]
gdp_dg_1=[]
dp_dt=[]
tdp_dt=[]
t_1_2=[]
t_1=[]
P_ISIP=[]
Flow=[]
Eff_Prs=[]
Term1=[]
Stiff=[]
h_fun=[]
Matrix=[]
Row=[]
log_dt=[]
log_tdpdt=[]
Input=[]

import csv, math
import numpy as np
import matplotlib.pyplot as plt

Reading()                           #Reading raw data and input parameters from SCV files

Young =Input[0]
Poisson=Input[1]
Poro =Input[2]
Comp=Input[3]
Fluid_Com=Input[4]
Fluid_Vis=Input[5]
Res_Pres_Inc=Input[6]
Frac_h=Input[7]
WSC=Input[8]
Inj_Total=Input[9]

Extract()                           #Extracting data

Compliance()                        #Calculating of G-function a and its dependant variables

CntPnt=Tangent()                           #Finding Shmin via tangent method
print(CntPnt,G_time[CntPnt-1],Prs_Shut[CntPnt-1], ShminT, gdp_dg_1[CntPnt])

h_function()

Results()

Eport_CSV()

fig,ax1=plt.subplots()
ax2=ax1.twinx()
ax3=ax1.twinx()
ax1.plot(G_time,Prs_Shut,'b')
ax2.plot(G_time[1:-1],dp_dg_1[1:],'r')
ax3.plot(G_time[1:-1],gdp_dg_1[1:],'g')
ax3.plot([0,G_time[CntPnt-1]],[0,gdp_dg_1[CntPnt]],'y')
ax3.plot([G_time[CntPnt-1],G_time[CntPnt-1]],[0,gdp_dg_1[CntPnt]],'y')
ax1.plot([0,G_time[CntPnt-1]],[ShminT,ShminT],'m')
ax1.plot([0,G_min_dpdg],[ShminC,ShminC],'c')
ax1.plot([G_min_dpdg,G_min_dpdg],[0,ShminC],'c')
ax1.set_xlim([0,60])
ax1.set_ylim([0, 10000])
ax2.set_ylim([0, 40])
ax3.set_ylim([0, 1100])

plt.show()