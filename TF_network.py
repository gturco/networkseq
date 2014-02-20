import commands
class Network(object):
    
    def __init__(self,TF,TF_id):
        self.children = []
        self.name = TF
        self.TF_id = TF_id

        
    def get_count_data(self,count_data):
        #expression_dic["AT1"]['Est.2364']["AT1G01050.1"]
        counts = [count_data[self.name]]
        for child in self.children:
            counts.append(count_data[child])
        return counts

def parse_network(TF_file):
    network_dict = {}
    for line in open(TF_file):
        ## use TF number here
        if line[0] == "": continue
        args = line.strip().split("\t")
        TF = args[3]
        line_id = args[2].replace("+AC0","")
        network_dict[line_id] = Network(TF,line_id)
        for target in args[3:]:
            if len(target) == 0: continue
            network_dict[line_id].children.append(target)
    return network_dict


def parse_countdata(expression_file,network_dict):
    ## messy dict lines = {"ATG1":{"1237":{"a":2300,"b":400,"d":0,"c":50}}}
    expression_dic = {}
    ex= open(expression_file)
    ex_data = ex.read().split("\n")
    genes = ex_data[0].strip().split(",")[1:]
    for line in ex_data[1:-1]:
        count_data = line.strip().split(",")
        TF_id = count_data[0].split(".")[1]
        TF_line_AT = network_dict[TF_id].name
        try:
            expression_dic[TF_line_AT][TF_id] = {}
        except KeyError:
            expression_dic[TF_line_AT] = {} 
            expression_dic[TF_line_AT][TF_id] = {}
            
        for i,gene in enumerate(genes):
            expression_dic[TF_line_AT][TF_id][gene.split(".")[0]] = count_data[i+1]

    #print expression_dic["AT1"]['Est.2364']["AT1G01050.1"]
    return expression_dic
            

def create_R_plot(TF, counts,WT):
    TF_c = counts[0]
    WTF_c = WT[TF.name]
    print TF.TF_id, float(TF_c)/float(WTF_c)
    for i,child in enumerate(TF.children):
        c = counts[i+1]
        wt_c = WT[child]
        out = open("tmp.txt",'wb')
        line = "{4}_{5},{6},Type\n{0},{1},WT\n{2},{3},Iducible\n".format(WTF_c,wt_c,TF_c,c,TF.name,TF.TF_id,child)
        out.write(line)
        out.close()
        x = commands.getstatusoutput("R CMD BATCH TF_plot.R")
        ##### num RNA CMD HERE

    
def main(TF_file, expression_file):
    
    network_dict = parse_network(TF_file)
    line_expression_data = parse_countdata(expression_file,network_dict)
    for key in network_dict.keys():
        TF = network_dict[key]
        print key
        ## MAKE Folder for this TF
        for ex_line in line_expression_data[TF.name].keys():
            print ex_line
            ex_line = line_expression_data[TF.name][ex_line]
            if TF.name == "WT" : continue
            counts = TF.get_count_data(ex_line)
            ### make folder in folder
            WT = line_expression_data["WT"]["Col-71"]
            create_R_plot(TF,counts,WT)
            print TF.name, TF.TF_id

main( "/Users/gturco/Desktop/code/Brady/md/LIBRARY INFO_RNA_SEQ_ESTRADIOL_01 3.csv","/Users/gturco/Desktop/code/Brady/md/all.txt")


