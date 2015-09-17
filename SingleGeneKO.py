############# Iterates through all genes in model to determine FBA outputs when each is knocked out. ##################

import os

#List of GPRs without reaction names from which list of genes is made
input1 = open('/gpfs/home/lmb5780/work/WM1788_N2fixing/K12_GPRS_NoRxnNames.txt','rU')
alllines1=input1.readlines()

## Make a list of genes (will iterate through)
geneList = []


for line in alllines1:
    #if line is empty, skip
    if line:
        if 'and' in line:
            #Processing for if line contains exclusively 'ands' or a combination of 'ands' and 'ors'
            #creates temporary lists of genes in individual GPRs
            tempGPRlist1 = line.split('and')
            for gene in tempGPRlist1:
                gene = gene.strip()
                if 'or' in gene:
                    #creates a new list for genes split by 'ands' which still contain 'or'. Removes from original list.
                    #Proceeds with stripping to just gene name
                    tempGPRlist2 = gene.split('or')
                    for gene in tempGPRlist2:
                        gene = gene.strip()
                        gene = gene.strip('()')
                        if gene not in geneList:
                            #Now that processing complete, add to list of genes if it has not already been added.
                            geneList.append(gene)
                #now moving on to case where line contains just "and"
                gene = gene.strip('()')
                if gene not in geneList:
                    #Now that split and strip complete, add to list of genes if it has not already been added.
                    geneList.append(gene)
        elif 'or' in line:
            #Processing for if line contains exclusively 'or'
            tempGPRlist3 = line.split('or')
            for gene in tempGPRlist3:
                gene = gene.strip()
                gene = gene.strip('()')
                if gene not in geneList:
                    #Now that processing complete, add to list of genes if it has not already been added.
                    geneList.append(gene)
        else:
            #GPR contains neither 'and' nor 'or' and is not empty.
            line = line.strip()
            line = line.strip('()')
            if line not in geneList:
             #Now that processing complete, add to list of genes if it has not already been added.
                geneList.append(line)

geneList.remove("")

#Testing to make sure the geneList was created correctly.
output1 = open('/gpfs/home/lmb5780/work/WM1788_N2fixing/ListOfGenes.txt', 'w')
for gene in geneList:
    output1.write(gene + '\n')
    
    
##Making a dictionary of reactions and their associated GPRs

#List of Reaction names with associated GPRS which will be used to make the reaction-GPR dictionary
input2 = open('/gpfs/home/lmb5780/work/WM1788_N2fixing/WM1788_GPR.txt','rU')
alllines2 = input2.readlines()

#Defining the reaction-GPR dictionary (and its copy for iterations later)
Reaction_GPR = dict()
Reaction_GPR_COPY = dict()

for line in alllines2:
    line = line.strip()
    rxn,GPR = line.split('\t')
    Reaction_GPR[rxn] = GPR
    
    
#for each gene
for gene in geneList:
    #list of rxns that will be False 
    knockoutList = []
    #creating a copy of Reaction_GPR dictionary that will be overwritten for each iteration/evaluation
    Reaction_GPR_COPY = Reaction_GPR.copy()
    #for each GPR
    for rxn in Reaction_GPR_COPY.keys():
        GPR = Reaction_GPR_COPY[rxn]
        #if gene present in GPR replace with "False" otherwise replace with "True"
        if GPR == "#N/A":
            Reaction_GPR_COPY[rxn] = 'True'
        if gene not in GPR:
           Reaction_GPR_COPY[rxn] = 'True'
        if gene in GPR:
            Reaction_GPR_COPY[rxn] = Reaction_GPR_COPY[rxn].replace(gene, 'False')
            for gene1 in geneList:
                if gene1 != gene and gene1 in GPR:
                    Reaction_GPR_COPY[rxn] = Reaction_GPR_COPY[rxn].replace(gene1, 'True')
        # eval
        # add to a dictionary of reactions with either false or true
        Reaction_GPR_COPY[rxn] = eval(Reaction_GPR_COPY[rxn])
        #Writing the list of knockouts
        if Reaction_GPR_COPY[rxn] == False:
            knockoutList.append(rxn)
    #Write knockout file for each gene
    geneOutput = open('/gpfs/home/lmb5780/work/WM1788_N2fixing/Knockouts_' + gene + '.txt', 'w')
    for reaction in knockoutList:
        geneOutput.write(reaction + '\n')
    #run GAMS to get max biomass, max and min ATP, max and min flavodoxin. 
    command = 'inputgms = open("/gpfs/home/lmb5780/work/WM1788_N2fixing/Knockouts.gms","r")'
    exec(command)
    full = inputgms.read()
    full = full.replace('$include knockouts','$include Knockouts_' + gene + '.txt')
    command = 'outputgms0=open("/gpfs/home/lmb5780/work/WM1788_N2fixing/max_biomass_WM1788_n2fix_KO.gms","w")'
    exec(command)
    outputgms0.write(full)
    outputgms0.close()
    cmd="gams max_biomass_WM1788_n2fix_KO.gms"
    os.system(cmd)
    #reading GAMS output
    input3 = open('/gpfs/home/lmb5780/work/WM1788_N2fixing/autotroph_max_report_KO.txt','rU')
    alllines3 = input3.readlines()
    for line in alllines3:
        if "Biomass" in line:
            biomass = line.split(':         ')[1]
        if "Minimum ATP" in line:
            min_ATP = line.split(':         ')[1]
        if "Maximum ATP" in line:
            max_ATP = line.split(':         ')[1]
        if "Minimum flavodoxin" in line:
            min_flavodoxin = line.split(':         ')[1]
        if "Maximum flavodoxin" in line:
            max_flavodoxin = line.split(':         ')[1]
    #add to file of genes knocked out with listed effects (max biomass, ATP, flavodoxin)
    output2 = open('/gpfs/home/lmb5780/work/WM1788_N2fixing/Gene_KO_EFFECTS','w')
    output2.write(gene + '\t' + 'Max Biomass' + '\t' + biomass + '\t' +  'Min ATP' + '\t' + min_ATP + '\t' + 'Max ATP' + '\t' + max_ATP + '\t' + 'Min Flavodoxin' + '\t' + min_flavodoxin + '\t' + 'Max Flavodoxin' + '\t' + max_flavodoxin + '\n')
    
    
    