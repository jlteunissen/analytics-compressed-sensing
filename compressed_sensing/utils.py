def generate_structures (RS_structures, ZB_structures):
    
    for i in range (len(RS_structures)):
        bulk_atoms = RS_structures[i].repeat(rep=3).get_chemical_symbols()
        bulk_positions = RS_structures[i].repeat(rep=3).get_positions()
        formula = RS_structures[i].get_chemical_formula()
        file = open("data/compressed_sensing/structures/RS_structures/"+formula+".xyz","w") 
        file.write ("%d\n\n"%54)
        for j in range (len(bulk_positions)):
            file.write (bulk_atoms[j])
            file.write ("\t%f\t%f\t%f\n"%(bulk_positions[j][0],bulk_positions[j][1],bulk_positions[j][2]))
        file.close()
        
    for i in range (len(ZB_structures)):
        bulk_atoms = ZB_structures[i].repeat(rep=3).get_chemical_symbols()
        bulk_positions = ZB_structures[i].repeat(rep=3).get_positions()
        formula = ZB_structures[i].get_chemical_formula()
        file = open("data/compressed_sensing/structures/ZB_structures/"+formula+".xyz","w") 
        file.write ("%d\n\n"%54)
        for j in range (len(bulk_positions)):
            file.write (bulk_atoms[j])
            file.write ("\t%f\t%f\t%f\n"%(bulk_positions[j][0],bulk_positions[j][1],bulk_positions[j][2]))
        file.close()
