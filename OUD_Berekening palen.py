# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import math
import numpy as np
import os
import matplotlib.pyplot as plt

#File input
paaljuk_folder = 'N02'
file = 'N02 - v1.4 - v2.1_incl_no_red CC2_inmeting.xlsx'
staaf_kolomnaam = 'Staaf'
aantal_palen = 6
paalnamen = False

# kies welke palen beschouwd dienen te worden. Hier dient de eerste paal uit de reeks genoemd te worden.
filter_toepassen = False
filter_staaf_output = ['S18', 'S68', 'S108', 'S113']
varianten_benaming = ['Ontgraven 1D', 'Ontgraven - JG10x', 'Ontgraven_JG10x_2D', 'Ontgraven_JG10x_hinge']

plot_legenda = True

buikmoment_in_staaf = False
beta_is_1 = ['S113'', S114', 'S115' , 'S116', 'S117']


#Input parameters voor berekening
t_dict={
   'Z14' : 12.5,
   'Z11' : 12.5,
   'Z08' : 12.5,
   'Z02' : 12.5,
   'Z05' : 12.5,
   'N11' : 12.5,
   'N08' : 12.5,
   'N08_rapportage' : 12.5,
   'N02' : 12.5,
   'N05' : 12.5,
   'Z15' : 8,
   'Z04' : 8,
   'Z06' : 8,
   'N02_RIUL' : 12.5,
   'Z06_2D' : 8,
   'Z04_2D' : 8,
   'Z04_JG' : 8,
   'Z07' : 8,
   'Z06_JG' : 8,
   'Z10' : 8
   } #Wanddikte
paaljuk = paaljuk_folder #naam palen
LengteInPaaljuk = 0 #Limitatie horizontale as
kniklengtes_dict = {
    'Z14' : {'S11':7000, 'S12':8000, 'S13':8500, 'S14':8500, 'S15':8500, 'S16':7500, 'S17':12500, 'S18':13500, 'S19':14000, 'S20':13500, 'S21':14500, 'S22':12500},
    'Z11' : {'S11':7000, 'S12':8000, 'S13':8000, 'S14':8500, 'S15':8000, 'S16':7000, 'S17':12500, 'S18':13000, 'S19':13000, 'S20':13500, 'S21':13000, 'S22':12000},
    'Z08' : {'S11':13000, 'S12':15000, 'S13':13000, 'S14':13500, 'S15':16000, 'S16':12500, 'S17':8000, 'S18':8000, 'S19':8000, 'S20':9000, 'S21':9000, 'S22':7500},
    'Z02' : {'S11':13500, 'S12':14500, 'S13':12000, 'S14':14000, 'S15':15000, 'S16':12500, 'S17':8000, 'S18':8000, 'S19':7500, 'S20':8500, 'S21':8500, 'S22':7500},
    'Z04' : {'S18':17000, 'S19':18000, 'S20':17500, 'S21':15000, 'S22':18000, 'S26':9000, 'S27':9500, 'S28':10000, 'S29':8000, 'S30':10000},
    'Z05' : {'S11':7500, 'S12':7500, 'S13':8500, 'S14':9000, 'S15':8500, 'S16':7000, 'S17':6500, 'S18':7000, 'S19':8000, 'S20':7500, 'S27':7500, 'S28':8000, 'S33':11500, 'S34':13500, 'S35':13500, 'S36':14000, 'S37':14500, 'S38':11000, 'S39':11000, 'S40':12500, 'S41':13500, 'S42':13000, 'S45':11500, 'S46':12500},
    'N11' : {'S11':6883, 'S12':6324, 'S13':6424, 'S15':7064, 'S16':5994, 'S17':7512, 'S18':13610, 'S19':12502, 'S20':12366, 'S21':12420, 'S22':11873, 'S23':13211, 'S28':12553, 'S29':11666, 'S30':11605, 'S31':11673, 'S32':11020, 'S33':12408},
    'N08' : {'S11':8500, 'S12':9000, 'S13':9000, 'S15':8500, 'S16':9000, 'S17':9000, 'S18':14000, 'S19':15500, 'S20':14000, 'S21':14500, 'S22':14500, 'S23':14500, 'S24':13500, 'S25':15000, 'S26':14000, 'S27':14000, 'S28':14500, 'S29':14000},
    'N02' : {'S11':8000, 'S12':9000, 'S13':8500, 'S14':8500, 'S15':8500, 'S16':8000, 'S17':14000, 'S18':17500, 'S19':14500, 'S20':14000, 'S21':14500, 'S22':17000, 'S23':13500, 'S24':17000, 'S25':14500, 'S26':13500, 'S27':14000, 'S28':16500, 'S34':12500, 'S35':14000, 'S36':13000, 'S37':14500, 'S38':16500, 'S39':15500},
    'N05' : {'S47':6961, 'S48':7850, 'S49':7709, 'S50':7728, 'S51':7843, 'S52':8916, 'S53':10367, 'S54':7863, 'S55':8060, 'S56':7500, 'S57':7175, 'S58':7156, 'S59':7933, 'S60':7974, 'S71':12249, 'S72':13087, 'S73':13054, 'S74':13342, 'S75':13836, 'S76':14885, 'S77':17170, 'S78':13064, 'S79':13893, 'S80':12689, 'S81':11890, 'S82':11645, 'S83':14518, 'S84':13687, 'S95':11759, 'S96':12595, 'S97':12554, 'S98':12821, 'S99':13283, 'S100':14297, 'S101':16535, 'S102':12583, 'S103':13362, 'S104':12220, 'S105':11465, 'S106':11242, 'S107':13868, 'S108':13180},
    'Z15' : {'P75-2':9600, 'P76-2':9600, 'P77-2':9600, 'P78-2':9600, 'P75-13':9600, 'P76-13':9600, 'P77-13':9600, 'P78-13':9600, 'P75-3':14500, 'P76-3':14500, 'P77-3':14500, 'P78-3':14500, 'P75-4':14600, 'P76-4':14600, 'P77-4':14600, 'P78-4':14600, 'P75-11':14600, 'P76-11':14600, 'P77-11':14600, 'P78-11':14600, 'P75-12':18600, 'P76-12':18600, 'P77-12':15900, 'P78-12':18600},
    'N02_RIUL' : {'S11':14500, 'S12':15500, 'S13':14000, 'S14':12500, 'S15':12500, 'S16':15500, 'S17':12000, 'S18':13500, 'S19':12500, 'S20':13500, 'S21':13500, 'S22':21000},
    'Z06' : {'S18':15576, 'S19':16517, 'S20':16370, 'S21':13782, 'S22':16904, 'S26':8721, 'S27':9322, 'S28':9674, 'S29':7914, 'S30':9827},
    'Z06_2D' : {'S18':15100, 'S19':16004, 'S20':15864, 'S21':13360, 'S22':16382, 'S26':16480, 'S27':17905, 'S28':17933, 'S29':15746, 'S30':17550, 'S49':10660, 'S50':11560, 'S51':12020, 'S52':10140, 'S53':11581, 'S63':15323, 'S64':16802, 'S65':16751, 'S66':14956, 'S67':16493},
    'Z04_2D' : {'S18':15140, 'S19':16030, 'S20':15813, 'S21':15140, 'S22':16360, 'S26':16324, 'S27':17700, 'S28':17705, 'S29':15555, 'S30':17293, 'S49':11762, 'S50':12750, 'S51':13100, 'S52':11167, 'S53':12710, 'S63':14598, 'S64':15983, 'S65':15853, 'S66':14205, 'S67':15641, 'S68':10880, 'S69':11595, 'S70':11847, 'S71':9790, 'S72':12085},
    'Z04_JG' : {'S68':10500, 'S69':11500, 'S70':11500, 'S71':10000, 'S72':11500},
    'Z07' : {'S8':8000, 'S9':7500, 'S10':7000, 'S11':7000, 'S12':14000, 'S13':13500, 'S14':13500, 'S15':14000, 'S16':11000, 'S17':10000, 'S18':10000, 'S19':12000, 'S20':11000, 'S21':10000, 'S22':10500, 'S23':11500},
    'Z06_JG' : {'S68':10500, 'S69':11500, 'S70':11500, 'S71':10000, 'S72':11500},
    'N08_rapportage' : {'S24':13500, 'S25':15000, 'S26':14000, 'S27':14000, 'S28':14500, 'S29':14000, 'S30':13500, 'S31':14500, 'S32':13500, 'S33':14000, 'S34':14000, 'S35':14000},
    'Z10' : {'S16':22000, 'S17':11500, 'S18':13000, 'S19':12000, 'S28':10000, 'S29':6500, 'S30':8500, 'S31':7500}
    }

geschoorde_palen = {
    'Z15' : ['S31', 'S32', 'S33', 'S34', 'S62', 'S63', 'S64', 'S65']
        }

#Kniklengtes uit SCIA
titel = file

kniklengtes = kniklengtes_dict[paaljuk_folder]
geschoorde_palen_juk = geschoorde_palen[paaljuk_folder] if paaljuk_folder in geschoorde_palen else []
t = t_dict[paaljuk_folder]

staven_filter_lst = []
for staaf_naam in filter_staaf_output:
    nummer = int(staaf_naam.replace('S',''))
    variant_lst = list(range(nummer, nummer + aantal_palen))
    variant_lst = ['S' + str(item) for item in variant_lst]
    staven_filter_lst.extend(variant_lst)


#Vastleggen van enkele vaste waardes voor de berekening. Deze waardes zijn overgenomen uit de berekening in het ontwerp van 2016, tenzij anders vermeld in de notitie
if t==8:
    EI=147337
    N_plRd=8271
    V_plRd=1257
    Alpha=0.21
    EA_ratio=2908/5842
    N_plrdbuis=4154
    beta=0.66
    M_plRd=674 # bij fy = 300
if t==12.5:
    EI=203552
    N_plRd=11629
    V_plRd=2800
    Alpha=0.21
    EA_ratio=4507/7361
    N_plrdbuis=7619 # was gelijk aan 6438, maar dat is berekend met fy=300 (incorrect)
    beta=0.66
    M_plRd=1227 # bij fy = 355

#Parameters voor het opslaan van de toetswaarde als figuur
plt.rcParams["figure.dpi"]=120
plt.rcParams["figure.figsize"]=[7.50*2,3.50*2]


#Opmaken van enkele lege sets en inlezen van de resultaten uit SCIA (mbhv resultaten.csv bestand). Toevoegen van enkele lege kolommen aan het dataframe
Staafnamen=set()
results=pd.read_excel(os.path.join(paaljuk_folder, file))
results["k"]=np.nan
results["n"]=np.nan
results["UC_N"]=np.nan
results["UC_V"]=np.nan
results["UC_NM"]=np.nan
results["L_cr"]=np.nan
results = results.rename(columns = {staaf_kolomnaam : 'Staaf'})


#Vastleggen van de locatie waar het opgeslagen moet worden
maploc=os.path.join(paaljuk_folder, f'{titel}_output')

if not os.path.exists(maploc):
        os.makedirs(maploc)
        
if filter_toepassen:
    results = results[results['Staaf'].isin(staven_filter_lst)]
    print(f'Resultaten gefilterd tot {staven_filter_lst}' )


#Elke rij van het dataframe gaan overlopen en de berekening gaan uitvoeren.
#De rijen "dx" bevatten vaak nog een "+" of "-" teken om aan te duiden of de snede aan de bovenzijde of onderzijde van het 1D-element is. Deze dienen eerst verwijderd te worden.
for index,row in results.iterrows():
    staaf = row['Staaf']
    if staaf not in kniklengtes:
        print(f'{row["Staaf"]} not found in kniklengtes')
        results = results.drop(index)
        continue
    
    if filter_toepassen and staaf not in staven_filter_lst:
        print(f'{row["Staaf"]} not found in staven filter')
        results = results.drop(index)
        continue

    if type(row['dx [mm]'])==str:

        if row['dx [mm]'][-1]=="+" :
            results.at[index,'dx [mm]']=float(row['dx [mm]'][:-1].replace(',','.'))
        elif row['dx [mm]'][-1]=="-" :
            results.at[index,'dx [mm]']=float(row['dx [mm]'][:-1].replace(',','.'))
        else:
            try:
                results.at[index,'dx [mm]']=float(row['dx [mm]'].replace(',','.'))
            except:
                results.at[index,'dx [mm]']=float(row['dx [mm]'])

    #Inlezen van kniklengtes, Xhi en Psi waarde berekenen per snede
    kniklengte=kniklengtes[staaf]
    N_cr=math.pi**2*EI/(kniklengte**2)*10**6
    Lambda=math.sqrt(N_plRd/N_cr)
    Psi=0.5*(1+Alpha*(Lambda**2-0.2)+Lambda**2)
    Xhi=1/(Psi+math.sqrt(Psi**2-Lambda**2))
    results.at[index,"L_cr"]=kniklengte
    if Xhi >1:
        Xhi=1
    if staaf in beta_is_1 and buikmoment_in_staaf:
        beta_calc = 1
    else:
        beta_calc = beta
    k=beta_calc/(1-abs(row['N [kN]'])/N_cr)
    if k<1:
        k=1
        results.at[index,'k']=k
    if k>1:
        results.at[index,'k']=k
    n=abs(row['N [kN]'])*EA_ratio/N_plrdbuis
    results.at[index,"n"]=n
    
    geschoord = True if staaf in geschoorde_palen_juk else False

    #Uitrekenen van de relevante UC
    results.at[index,"Xhi"] = Xhi
    results.at[index,"UC_N"]=abs(row["N [kN]"])/(Xhi*N_plRd)
    results.at[index,"N [kN]"]=abs(row["N [kN]"])
    V_tot=math.sqrt(row["Vy [kN]"]**2+row["Vz [kN]"]**2)
    results.at[index,"VEd [kN]"] = V_tot
    results.at[index,"UC_V"]=V_tot/V_plRd
    Med = math.sqrt(row["My [kNm]"]**2+row["Mz [kNm]"]**2)
    if geschoord:
        results.at[index,"MEd [kNm]"] = Med
        results.at[index,"MEd_k [kNm]"] = Med
        results.at[index,"UC_NM"]=Med/M_plRd
    else:
        M_tot=Med*k
        results.at[index,"MEd [kNm]"] = Med
        results.at[index,"MEd_k [kNm]"] = M_tot
        M_NyRd=M_plRd*(1-n**(1.7))
        results.at[index,"UC_NM"]=M_tot/M_NyRd
    Staafnamen.add(staaf)

pivot_for_uc_table = pd.pivot_table(results, values=["N [kN]", 'UC_N', 'VEd [kN]', 'UC_V', 'MEd [kNm]', 'MEd_k [kNm]', 'UC_NM'], index='Staaf', aggfunc=np.max).round(2)
# add thickness column
pivot_for_uc_table['Buisdikte [mm]'] = t
# add kniklengtes
df = pd.DataFrame.from_dict(kniklengtes, orient='index', columns = ['Kniklengte [m]']).apply(lambda x: x/1000).round(1)
pivot_for_uc_table= pivot_for_uc_table.join(df)


pivot_for_uc_table = pivot_for_uc_table[['Buisdikte [mm]', 'Kniklengte [m]', 'N [kN]', 'MEd [kNm]', 'MEd_k [kNm]', 'VEd [kN]',  'UC_N', 'UC_NM', 'UC_V']]
pivot_for_uc_table[['N [kN]', 'MEd [kNm]', 'MEd_k [kNm]', 'VEd [kN]']] = pivot_for_uc_table[['N [kN]', 'MEd [kNm]', 'MEd_k [kNm]', 'VEd [kN]']].round(0).astype(int)

# sort according to dict order
pivot_for_uc_table = pivot_for_uc_table.loc[results.Staaf.unique()]

#Opslaan van het bestand naar een xlsx bestand zodat de output herleidbaar is
excel_dir = maploc +"/"+titel +" - output.xlsx"

#Create or overw empty excel file to fill
pd.DataFrame().to_excel(excel_dir)


with pd.ExcelWriter(excel_dir, mode='a', if_sheet_exists='replace') as writer:  
    results.to_excel(writer, sheet_name = 'data_tabel', index=False)
    pivot_for_uc_table.to_excel(writer, sheet_name='UC_tabel',index=True)


#sort staafnamen
staaflst = list(Staafnamen)
staaflst.sort()

if paalnamen:
    variant_nummers = list(set([paalnaam.split('-')[1] for paalnaam in staaflst]))
    variant_nummers.sort()
    staaflst0 = []
    for i in variant_nummers:
        staaflst1 = [staaf for staaf in staaflst if f'-{i}' in staaf]
        staaflst0.extend(staaflst1)
    staaflst = staaflst0
else:
    #fix sort with different lengths
    staaflst0 = [staaf for staaf in staaflst if len(staaf) == 2]
    staaflst1 = [staaf for staaf in staaflst if len(staaf) == 3]
    staaflst2 = [staaf for staaf in staaflst if len(staaf) == 4]

    staaflst0.extend(staaflst1)
    staaflst0.extend(staaflst2)
    staaflst = staaflst0

#Opmaken van plots
aantal_plots = min(len(kniklengtes.keys()), len(staaflst))
aantal_varianten = len(staven_filter_lst) / aantal_palen if filter_toepassen else aantal_plots / aantal_palen
aantal_varianten = max(int(aantal_varianten), 2)
fig, axs = plt.subplots(aantal_varianten, aantal_palen)

for i, staaf in enumerate(Staafnamen):
    y_b=[]
    y=results[results.Staaf==staaf]
    locations=set()
    for index,row in y.iterrows():
        locations.add(row["dx [mm]"])

    #Vastleggen van x-as
    for location in locations:
        y_a=y[y["dx [mm]"]==location]
        y_b.append([float(location),max(y_a.UC_N),max(y_a.UC_V),max(y_a.UC_NM)])

    locations=list(locations)
    for index in range(0,len(locations)):
        locations[index]=float(locations[index])
    locations=sorted(locations)

    #Toevoegen van datase
    y_tot=pd.DataFrame(y_b,columns=["x","UC_N","UC_V","UC_NM"])
    y_tot=y_tot.sort_values("x")

    #Opmaken van grafiek
    paalnr = staaflst.index(staaf) + 1
    rij = int(np.ceil(paalnr/aantal_palen) - 1)
    kolom = int(paalnr - rij*aantal_palen) - 1
    
    print(staaf, paalnr, rij, kolom)
    
    uc_max = y[['UC_V', 'UC_N', 'UC_NM']].max().max().round(2)
    
    plot_title = f'{staaf} - {str(uc_max)}'
    y_tot.plot(x="x", ax=axs[rij, kolom], title = plot_title, legend = False)
    axs[rij, kolom].set_xlim(0,locations[-1]-LengteInPaaljuk)
    axs[rij, kolom].set_ylim(0,1.5)
    axs[rij, kolom].set_xlabel("Positie [mm]")
    axs[rij, kolom].set_ylabel("UC [-]")
    axs[rij, kolom].hlines(y=1,xmin=0,xmax=locations[-1],color='r', linestyle='dashed')
    
    if filter_toepassen and kolom == 0:
        variant_benaming = varianten_benaming[filter_staaf_output.index(staaf)]
        axs[rij, kolom].set_title(f'{variant_benaming} \n {staaf} - {str(uc_max)}')

for ax in axs.flat:
    ax.label_outer()

handles, labels = ax.get_legend_handles_labels()
# fig.figlegend()
if plot_legenda:
    fig.legend(handles, labels, loc='center right')
fig.suptitle(titel)
plt.subplots_adjust(hspace = 0.3)

fig.savefig(maploc+"/" + os.path.splitext(file)[0] +".pdf")


#Opstellen van tekstbestand met rekenparameters
Rekenparameters=maploc +"/"+paaljuk +"rekenparameters.txt"
parameterfile=open(Rekenparameters,'w')
parameterfile.write(titel)
parameterfile.write(" \n")
for element in kniklengtes:
    parameterfile.write("Kniklengte staaf \t"+element+"\t \t"+str(kniklengtes[element])+" mm")
    parameterfile.write("\n")
parameterfile.write("Buisdikte [mm]"+"\t"+str(t)+"\n")
parameterfile.write("Lengte in paaljuk [mm]"+"\t"+str(LengteInPaaljuk)+"\n")
parameterfile.write("EI [10^9 Nmm²]"+"\t"+str(EI)+"\n")
parameterfile.write("N_plRd [kN]"+"\t"+str(N_plRd)+"\n")
parameterfile.write("V_plRd [kN]"+"\t"+str(V_plRd)+"\n")
parameterfile.write("EA_buis / EA_totaal [-]"+"\t"+str(EA_ratio)+"\n")
parameterfile.write("N_plRd,buis [kN]"+"\t"+str(N_plrdbuis)+"\n")
parameterfile.write("M_plRd [kNm]"+"\t"+str(M_plRd)+"\n")
parameterfile.write("N_cr [kN]"+"\t"+str(N_cr)+"\n")
parameterfile.write("Lambda [-]"+"\t"+str(Lambda)+"\n")
parameterfile.write("Xhi [-]"+"\t"+str(Xhi)+"\n")
parameterfile.write("\n \n")    
parameterfile.write(str(kniklengtes))

parameterfile.close()
