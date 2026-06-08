import os, linecache
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import MultipleLocator
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #Only needed if modifying defaults

np.set_printoptions(threshold=np.inf)
np.set_printoptions(suppress=True)
np.set_printoptions(linewidth=1000)#不自动换行
np.set_printoptions(precision=6) 

#config = {"mathtext.fontset":'stix'}
#rcParams.update(config)


def disorted_NAO_OVC(out_file,excited_file,mol_file):
    mol = Chem.MolFromPDBFile(mol_file)
    ssr = Chem.GetSymmSSSR(mol)

    f = open(out_file,'r')
    line = f.readlines()
    occ = []
    for index,lines in enumerate(line):
        if "Charge =  0 Multiplicity = 1" in lines:
            coord_start = index+1
        if "NQM=" in lines:
            natom_start = index
        if "basis functions" in lines:
            basis_start = index
        if "Alpha  occ. eigenvalues" in lines:
            occ_start = index
            occ.append(occ_start)
        if "NATURAL POPULATIONS:  Natural atomic orbital occupancies " in lines:
            NAO_start = index+4
        if "NAO_search" in lines:
            NAO_value_start = index+2
        if "Center     Atomic      Atomic" in lines:
            correct_coord = index+3       

    ### 
    natom = np.loadtxt(out_file,skiprows=natom_start,max_rows=1,usecols=3,dtype='int')
    nbasis = np.loadtxt(out_file,skiprows=basis_start,max_rows=1,usecols=0,dtype='int')
    atom_name = np.loadtxt(out_file,dtype='str',skiprows=coord_start,max_rows=natom,usecols=0)
    atom_coord = np.loadtxt(out_file,skiprows=correct_coord,max_rows=natom,usecols=(3,4,5))

    occ_last_line = np.loadtxt(out_file,skiprows=occ[-1],max_rows=1,dtype='str')
    HOMO = (len(occ)-1)*5+len(occ_last_line)-4
    LUMO = (len(occ)-1)*5+len(occ_last_line)-4+1


    NAO_lang=np.loadtxt(out_file,skiprows=NAO_start,max_rows=nbasis+natom-1,dtype='str',comments='#',usecols=3)
    NAO_type_1=np.loadtxt(out_file,skiprows=NAO_start,max_rows=nbasis+natom-1,dtype='str',comments='#',usecols=4)
    NAO_type = [i[0:-1] for i in NAO_type_1]


    #键连关系
    not_H_atom = []
    for i in range(len(atom_name)):
        if atom_name[i] != 'H':
            not_H_atom.append(i)            
    not_H_coord = atom_coord[not_H_atom]

    ring_atom = []
    all_ring = []
    for ring in ssr:
        all_ring.append(list(ring))
        for j in list(ring):
            ring_atom.append(j)
    ring_atom = list(set(ring_atom))
    #不在环上原子
    not_in_ring = []
    for i in range(len(not_H_atom)):
        if i not in ring_atom:
            not_in_ring.append(i)

    bond_compute=[]
    bond_show=[]
    for i in range(len(not_H_coord)):
        for j in range(len(not_H_coord)):
            if j > i:
                distance = np.sqrt(np.sum((np.array(not_H_coord[i]) - np.array(not_H_coord[j])) ** 2))
                if 1.0 < distance < 2.0:
                    bond_compute.append([i,j])
                    bond_show.append([not_H_atom[i]+1,not_H_atom[j]+1])

    ###孤悬在外的原子###
    single_atom=[]
    for i in not_in_ring:
        for j in np.array(bond_compute):
            if i in j:
                single_atom.append([i,j])
    single_atom_bondatom=[]
    for i in single_atom:
        for j in i[1]:
            if j != i[0]:
                single_atom_bondatom.append([i[0],j])
    single_atom_ringatom = []  #孤悬原子所连接的环 
    for i in single_atom_bondatom:
        for j in all_ring:
            if i[1] in j:
                single_atom_ringatom.append([i[0],j])
                
    #tBu原子!!!

    tBu_atom = []
    for i in not_in_ring:
        if i not in [x[0] for x in single_atom_ringatom]:
            tBu_atom.append(i)
            
    #环上共用的原子
    share_atom_show = []
    for i in np.array(not_H_atom)+1:
        share_atom_1 = []
        for j in bond_show:
            if i == j[0] or i == j[1]:
                share_atom_1.append(i)
        if len(share_atom_1) == 3:
              share_atom_show.append(i)
    #环上共用的原子                             
    share_atom_compute = []   
    share_atom_4 = []
    share_atom_4_acompany = []#
    for i in range(len(not_H_atom)):            
        share_atom_2 = [] 
        share_atom_2_acompany = []
        share_atom_2_acompany.append(i)
        for k in bond_compute:
            if i == k[0]:
                share_atom_2.append(i)
                share_atom_2_acompany.append(k[1])
            if i == k[1]:
                share_atom_2.append(i)
                share_atom_2_acompany.append(k[0])          

        if len(share_atom_2) == 3:
            share_atom_compute.append(i)
            share_atom_4_acompany.append(share_atom_2_acompany)


    #价层轨道
    coeff_num = []
    for i in range(nbasis):
        if NAO_lang[i] == 'px' and NAO_type[i] == 'Val':
            coeff_num.append(i)
        if NAO_lang[i] == 'py' and NAO_type[i] == 'Val':
            coeff_num.append(i)        
        if NAO_lang[i] == 'pz' and NAO_type[i] == 'Val':
            coeff_num.append(i)

    def nao_orb_coeff(orbital):        
        NAO_value_num = int((np.float64(orbital)-0.1)/8)*8+1
        blank_num = int(9)-int(len(str(NAO_value_num)))
        NAO_search = 'NAO'+' '*blank_num+str(NAO_value_num)   

        f = open(out_file,'r')
        line = f.readlines()
        for index,lines in enumerate(line):
            if str(NAO_search) in lines:
                NAO_value_start = index+2

    #原子轨道系数            
        NAO_value_orb = linecache.getlines(out_file)[NAO_value_start: NAO_value_start+nbasis]
        orb_value = []
        for line in NAO_value_orb:
            coef = line[17+8*int(orbital-NAO_value_num):17+8*int(orbital-NAO_value_num+1)]
            orb_value.append(coef)

    #价层系数
        orb_value_val = []
        for i in coeff_num:
            orb_value_val.append(np.float64(orb_value[i]))

        return np.array(orb_value_val)

    #校正的价层系数
    def adjuest_nao_orb_coeff(orbital):
        #环上的系数
        all_ring_coeff = []    
        for ring in ssr:
            ring_coord = not_H_coord[ring]
            center = ring_coord.mean(axis=0) 
            coordR = ring_coord - center #将分子移到中心
            u,sigma,v = np.linalg.svd(coordR)
            pnorm = v[2] #法向量
            if pnorm[2] < 0:
                pnorm=pnorm*-1
            for i in ring:
                all_ring_coeff.append([i,np.dot(pnorm,nao_orb_coeff(orbital).reshape(-1,3)[i])])

        del_list=[]
        for i in share_atom_4_acompany:
            del_list.append(i[0])

        new_ring_coeff = []
        for i in all_ring_coeff:    
            if i[0] not in del_list:
                new_ring_coeff.append(i)

        ####共享原子系数
        all_point_coeff = []
        for i in share_atom_4_acompany:
            single_4_coord = not_H_coord[i]
            center = single_4_coord.mean(axis=0) 
            coordR = single_4_coord - center #将分子移到中心
            u,sigma,v = np.linalg.svd(coordR)
            pnorm = v[2] #法向量 
            if pnorm[2] < 0:
                pnorm=pnorm*-1

            all_point_coeff.append([i[0],np.dot(pnorm,nao_orb_coeff(orbital).reshape(-1,3)[i[0]])])

        #孤悬原子系数
        single_atom_coeff = []
        for i in single_atom_ringatom:
            ringatom_coord = not_H_coord[i[1]]
            center = ringatom_coord.mean(axis=0) 
            coordR = ringatom_coord - center #将分子移到中心
            u,sigma,v = np.linalg.svd(coordR)
            pnorm = v[2] #法向量 
            if pnorm[2] < 0:
                pnorm=pnorm*-1
            single_atom_coeff.append([i[0],np.dot(pnorm,nao_orb_coeff(orbital).reshape(-1,3)[i[0]])])
            
        ###tBu原子系数1!1！1
        tBu_atom_coeff = []
        for i in tBu_atom:
            tBu_atom_coeff.append([i,np.dot([0,0,1],nao_orb_coeff(orbital).reshape(-1,3)[i])])

        #所有非H原子系数
        adjuested_coeff = new_ring_coeff+all_point_coeff+single_atom_coeff+tBu_atom_coeff
        adjuested_coeff_listed = sorted(adjuested_coeff,key=(lambda x:x[0]))

        adjuested_orb_value_val = []
        for i in adjuested_coeff_listed:
            adjuested_orb_value_val.append(i[1])
        return np.array(adjuested_orb_value_val)

    #print(adjuest_nao_orb_coeff(HOMO))
    #激发态跃迁信息
    e_f = open(excited_file,'r')
    e_line = e_f.readlines()
    e_start = []
    e_end = []

    for index,lines in enumerate(e_line):
        if "Excited State   1:" in lines:
            excited_start = index+1
            e_start.append(excited_start)
        if "This state for optimization and/or second-order correction." in lines:
            excited_end = index
            e_end.append(excited_end) 

    ###       
    trans_information = linecache.getlines(excited_file)[e_start[-1]:e_end[-1]]

    trans_orb_1 = []
    trans_coeff_1 = []
    for line in trans_information:
        init_orb = line[0:9]

        fina_orb = line[11:16]
        coeff_value_1 = line[22:35]
        trans_orb_1.append(init_orb)
        trans_orb_1.append(fina_orb)
        trans_coeff_1.append(coeff_value_1)

    trans_orb = np.array(trans_orb_1).reshape(-1,2)
    trans_coeff= np.array(np.float64(trans_coeff_1)).reshape(-1,1)


    #对应跃迁轨道的eta,phi系数
    def nao_coeff(inital,final):
        eta_inital_orb=[]
        eta_final_orb=[]
        phi_orb = []
        adjuest_coeff_out_inital = adjuest_nao_orb_coeff(int(inital))
        adjuest_coeff_out_final = adjuest_nao_orb_coeff(int(final))
        for i in bond_compute:
            eta_inital_orb.append((adjuest_coeff_out_inital[i[0]]*adjuest_coeff_out_inital[i[1]])**2)

            eta_final_orb.append((adjuest_coeff_out_final[i[0]]*adjuest_coeff_out_final[i[1]])**2)    
            phi_orb.append(adjuest_coeff_out_inital[i[0]]*adjuest_coeff_out_final[i[0]]*\
                           adjuest_coeff_out_inital[i[1]]*adjuest_coeff_out_final[i[1]])
        return sum(eta_inital_orb), sum(eta_final_orb), sum(phi_orb), np.array(eta_inital_orb), np.array(eta_final_orb), np.array(phi_orb)



    eta_inital_singlet_value=[]
    eta_final_singlet_value=[]
    phi_singlet_value=[]
    eta_value=[]
    phi_value=[]

    nao_coeff_out = []#跃迁系数只算一次
    for i in range(len(trans_orb)):
        nao_coeff_out.append(nao_coeff(trans_orb[i][0],trans_orb[i][1]))
    # eta,phi系数 乘 跃迁系数   
    for i in range(len(trans_orb)):
        eta_inital_value = nao_coeff_out[i][3]*(trans_coeff[i]*np.sqrt(2))**2
        eta_inital_singlet_value.append(eta_inital_value)

        eta_final_value = nao_coeff_out[i][4]*(trans_coeff[i]*np.sqrt(2))**2
        eta_final_singlet_value.append(eta_final_value)

        phi_value_1 = nao_coeff_out[i][5]*(trans_coeff[i]*np.sqrt(2))**2*-1
        phi_singlet_value.append(phi_value_1)

        eta_1= (nao_coeff_out[i][0]+nao_coeff_out[i][1])
        eta_2=eta_1*(trans_coeff[i]*np.sqrt(2))**2

        phi = nao_coeff_out[i][2]*(trans_coeff[i]*np.sqrt(2))**2*-1
        eta_value.append(eta_2)
        phi_value.append(phi)

    header = 'eta_value: '+str(np.float64(eta_value))+'\n'+'phi_value: '+str(np.float64(phi_value))+'\n'\
    +'Bond     ' +'eta_inital_value*N   '+'eta_final_value*N   '+'phi_value*N '
    print_data = []
    print_data.append(bond_show)
    for i in range(len(trans_orb)):  
        print_data.append(np.array(eta_inital_singlet_value)[i])#.reshape(len(trans_orb),-1))
    for i in range(len(trans_orb)):
        print_data.append(np.array(eta_final_singlet_value)[i])#.reshape(len(trans_orb),-1))
    for i in range(len(trans_orb)):    
        print_data.append(np.array(phi_singlet_value)[i])#.reshape(len(trans_orb),-1))
    

    with open(str(out_file[0:-4])+'.dat','ab') as f:
        np.savetxt(f,np.array(print_data).T,header=header,fmt=' '.join(['%-10s'] + ['%16.7f']*3*len(trans_orb))) 
        np.savetxt(f,np.hstack(((np.array(not_H_atom)+1).reshape(-1,1),\
np.array(adjuest_nao_orb_coeff(HOMO).reshape(-1,1)))),header='\n'+'\n'+'HOMO',fmt=' '.join(['%-10s'] +['%16.5f']))
        np.savetxt(f,np.hstack(((np.array(not_H_atom)+1).reshape(-1,1),\
np.array(adjuest_nao_orb_coeff(LUMO).reshape(-1,1)))),header='\n'+'\n'+'LUMO',fmt=' '.join(['%-10s'] +['%16.5f']))
    nao_orb_coeff
    f.close()
    
#os.chdir('F:/xmu/ARTADF/OL/PBE0/NAO-revise/basis')
#os.chdir('C:/Users/wanli/Desktop/reorganization/anthracene/benzene')
#out_file = '2d-s1-nao.log'
#excited_file = '2d-s1-vert.log'
#mol_file = '2d-s1-nao.mol'
#disorted_NAO_OVC('2d-s0-nao.log','2d-s1-vert.log','2d-s1-nao.mol')




logFiles = []
for filename in os.listdir('.'): 
    if filename.endswith('-nao.log'):
        logFiles.append(filename)
        
for i in logFiles:
    disorted_NAO_OVC(i,i[0:-10]+'s1-vert.log',i[0:-8]+'.pdb')