import sys

                                                               ######################### PREPARE ALL THE FUNCTIONS #####################
#1. FUNCTIONS FOR STATISTICS
  
stats = {"reads": 0, "bases": 0, "A": 0, "C": 0, "G": 0, "T": 0, "N": 0}
trimmed_stats = {"bases": 0, "A": 0, "C": 0, "G": 0, "T": 0, "N": 0}

#### STATS UPDATE

  #### stats update 
  
def stats_update(stats,sequence):
 stats["reads"] += 1 #sum one in each round
 stats["bases"] += len(sequence) #sum the lenght of the sequence in each round
 for nt in sequence:
  stats[nt] += 1 #sum 1 in each round to each of the bases
  
  
  #### Trimmed stats update 
 
def trimmed_stats_update(trimmed_stats,sequence,value):
  #sum the lenght of the sequence in each round
 
 if "--trim-right" in arguments and not "--trim-left" in arguments:
  right_modif_sequence= sequence[value:]
  trimmed_stats["bases"] += len(right_modif_sequence)
  for nt in right_modif_sequence:
   trimmed_stats[nt] += 1 
    
 elif "--trim-left" in arguments and not "--trim-right" in arguments:
  left_modif_sequence= sequence[:value]
  trimmed_stats["bases"] += len(left_modif_sequence)
  for nt in left_modif_sequence:
   trimmed_stats[nt] += 1 
 
def trimmed_stats_update_both(trimmed_stats,sequence,value_left,value_right):
  right_modif_sequence= sequence[value_right:]
  left_modif_sequence= sequence[:value_left]
  trimmed_stats["bases"] += (len(right_modif_sequence) + len(left_modif_sequence))
  for nt in left_modif_sequence:
   trimmed_stats[nt] += 1 
  for nt in right_modif_sequence:
   trimmed_stats[nt] += 1  
   
#### STATS SUMMARY  
  
def print_summary(arguments,stats,trimmed_stats=None): #you print a summary with the arguments and the stats and by default with a trimmed=none
# Translate operation
 operation = arguments["--operation"] #for the argument operation print the corresponding output
 if operation == "rc": 
  operation = "reversed-complemented"
 elif operation == "trim":
  operation = "hard-trimmed"
 else:
  operation = "processed " #if the user do´t put anything it will just process the files
  
# Print summary
 print("File '%s' has been successfully %s ('%s')" %
(arguments["--input"],operation,arguments["--output"]))
 print("Summary:");
 if arguments["--operation"]=="adaptor-removal" and argument=="--adaptor":
  print("\t%d reads processed" % stats["reads"])
  print("\t%d bases processed (%d%% A, %d%% C, %d%% G, %d%% T, %d%% N)" % (stats["bases"],(stats["A"]*100/stats["bases"]),(stats["C"]*100/stats["bases"]),(stats["G"]*100/stats["bases"]),(stats["T"]*100/stats["bases"]),(stats["N"]*100/stats["bases"])))
  print("\t%d adaptors-found" % num_adaptors)
 
 else:
  print("\t%d reads processed" % stats["reads"])
  print("\t%d bases processed (%d%% A, %d%% C, %d%% G, %d%% T, %d%% N)" % (stats["bases"],(stats["A"]*100/stats["bases"]),(stats["C"]*100/stats["bases"]),(stats["G"]*100/stats["bases"]),(stats["T"]*100/stats["bases"]),(stats["N"]*100/stats["bases"])))

 if trimmed_stats != None:
  print("\t%d bases trimmed (%d%% A, %d%% C, %d%% G, %d%% T, %d%% N)" % (trimmed_stats["bases"],(trimmed_stats["A"]*100/trimmed_stats["bases"]),(trimmed_stats["C"]*100/trimmed_stats["bases"]),
      (trimmed_stats["G"]*100/trimmed_stats["bases"]),(trimmed_stats["T"]*100/trimmed_stats["bases"]),(trimmed_stats["N"]*100/trimmed_stats["bases"])));



#2. FUNCTIONS TO MODIFY THE SEQUENCE

complement ={"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"} 


####REVERSE COMPLEMENT FUNCTION

def reverse_complement (sequence):
 modif_sequence="".join([complement[nt] for nt in sequence[::-1]])  
 return modif_sequence

#####FUNCTION TO REMOVE THE ADAPTOR

      #### Function to remove the adaptor from the sequence
      
def adaptor_removal(sequence,adaptor):
 count_adaptor=0
 if sequence[0:len(adaptor)]== adaptor:
  count_adaptor+=1
  begin = 0
  while begin < len(adaptor):
   begin += 1
  end = len(sequence) - 1
  return sequence[begin:end+1]
 else:
  return sequence
   

      #### Function to remove the adaptor from the qualities
      
def adaptor_removal_qualities(qualities,adaptor):
 return qualities[len(adaptor):] 
  

# 3. PARSING ARGUMENTS

arguments ={} #I create an empty dictionary which would be filled with the things indicated from the user 

i=1 
while i <len(sys.argv):
 argument =sys.argv [i]
 if argument != "--input" and argument !="--output" and argument !="--operation" and argument !="--adaptor" and argument !="--trim-right" and argument !="--trim-left":
   print("Unknown argument" "%s" % argument)
   exit(1)
 arguments[argument]=sys.argv[i+1]
 i +=2


# 4. FUNCTION TO OPEN THE INPUT FILE AND OUTPUT FILE 

### to open the input file

def open_input(file_name):
 input_file = open(file_name,"rt")
 line = input_file.readline()
 if line[0]==">":
  input_format = "FASTA"
 elif line[0]=="@":
  input_format = "FASTQ"
 else:
  print("File '%s' is neither a FASTQ/FASTA file" % file_name)
  exit(1)
 input_file.seek(0) 
 return (input_file,input_format) 
 
### to open the output file 

def open_output(file_name):
 output_file = open(file_name,"wt")
 return output_file
 
 
 

                                                             ######################### READ THE FILES #####################

(input_file, input_format)= open_input(arguments["--input"]) # the function open_input returns an input_file and an input format when you give the argument "input"
output_file = open_output(arguments["--output"]) # the function open_ouput returns an output file with the name the user indicates in the argument "output"
num_adaptors=0 #a count of the number of adaptors before the loop
while True: 
 tag =input_file.readline()[1::].rstrip("\n") #of the input file that you have opened read all the lines and separate all of them with a jump of line
 if not tag: break #if there are not lines stop the loop
 
 sequence=input_file.readline().rstrip("\n") #read the sequence, from line number 1
 
 if arguments["--operation"]=="rc":  #REVERSE COMPLEMENT:
  modif_sequence=reverse_complement(sequence)# generate a reverse complement sequence
    
 elif arguments["--operation"]=="adaptor-removal" and argument=="--adaptor":  #ADAPTOR REMOVAL:    
  adaptor=arguments["--adaptor"]
  adaptor=adaptor.upper()
  modif_sequence= adaptor_removal(sequence, adaptor) # remove the adaptor
  beg_seq=sequence[0:len(adaptor)] 
  num_adaptors += beg_seq.count(adaptor) #count only if there is an adaptor in the beggining of the sequence  
  
 elif arguments["--operation"]=="trim" and "--trim-right" in arguments and not "--trim-left" in arguments: #TRIM:
  value= len(sequence)-int(arguments["--trim-right"])
  modif_sequence= sequence[0:value] #we take the bases from the beggining of the sequence until the right position that has to be trimmed     
 elif arguments["--operation"]=="trim" and "--trim-left" and not "--trim-right" in arguments:
  value= int(arguments["--trim-left"])
  modif_sequence= sequence[value:] #we take the bases from the left position of the end of the trimming until the end of the sequence 
 elif arguments["--operation"]=="trim" and "--trim-left" and "--trim-right" in arguments:
  value_left= int(arguments["--trim-left"])
  value_right= len(sequence)- int(arguments["--trim-right"])
  modif_sequence= sequence[value_left:value_right]   
    
 if input_format=="FASTQ": #FORMAT: if input format is a fastq, according to what I have defined in open_input function
  input_file.readline() #read the +
  qualities = input_file.readline() #read qualities
  
  if arguments["--operation"]=="rc":
   output_file.write("@%s\n%s\n+%s\n" % (tag,modif_sequence,qualities[::-1]))  
  
  elif arguments["--operation"]=="adaptor-removal" and argument=="--adaptor":
   if sequence[0:len(adaptor)]== adaptor and len(adaptor)==len(sequence):
    output_file.write("")    
   else: 
    if sequence[0:len(adaptor)]== adaptor:
     qualities= adaptor_removal_qualities(qualities,adaptor)
    output_file.write("@%s\n%s\n+\n%s" % (tag,modif_sequence,qualities))
      
  elif arguments["--operation"]=="trim" and "--trim-right" and not "--trim-left" in arguments :  
   if int(arguments["--trim-right"])>=len(sequence):  
    output_file.write("")  
   else:
    qualities= qualities[0:value]
    output_file.write("@%s\n%s\n+\n%s\n" % (tag,modif_sequence,qualities))           
  elif arguments["--operation"]=="trim" and "--trim-left" and not "--trim-right" in arguments:
   if value>=len(sequence):  
    output_file.write("") 
   else:   
    qualities= qualities[value:]
    output_file.write("@%s\n%s\n+\n%s" % (tag,modif_sequence,qualities))
  elif arguments["--operation"]=="trim" and "--trim-left" and "--trim-right" in arguments:
   if (int(arguments["--trim-right"]) + value_left)>=len(sequence):  
    output_file.write("") 
   else:
    qualities= qualities[value_left:value_right]
    output_file.write("@%s\n%s\n+\n%s\n" % (tag,modif_sequence,qualities))  
  
  elif arguments["--operation"]=="trim" and not "--trim-left" and not "--trim-right" in arguments:
   print("Please, introduce a --trim-left or --trim-right value to trim") 
   output_file.write("@%s\n%s\n+\n%s" % (tag,sequence,qualities)) 
  
 if input_format=="FASTA":
  if arguments["--operation"]=="adaptor-removal" and argument=="--adaptor":
   if sequence[0:len(adaptor)]== adaptor and len(adaptor)==len(sequence):
    output_file.write("")   
   else:
    output_file.write(">%s\n%s\n" % (tag,modif_sequence))      
  elif arguments["--operation"]=="trim" and "--trim-right" and not "--trim-left" in arguments :  
   if int(arguments["--trim-right"])>=len(sequence):  
    output_file.write("")  
   else:
    output_file.write(">%s\n%s\n" % (tag,modif_sequence)) 
  elif arguments["--operation"]=="trim" and "--trim-left" and not "--trim-right" in arguments:
   if value>=len(sequence):  
    output_file.write("") 
   else:   
    output_file.write(">%s\n%s\n" % (tag,modif_sequence)) 
  elif arguments["--operation"]=="trim" and "--trim-left" and "--trim-right" in arguments:
   if (int(arguments["--trim-right"]) + value_left)>=len(sequence):  
    output_file.write("") 
   else:
    output_file.write(">%s\n%s\n" % (tag,modif_sequence)) 
  else:
   output_file.write(">%s\n%s\n" % (tag,modif_sequence)) 
     
# Update stats
 stats_update(stats,sequence)
 if arguments["--operation"]=="trim" and ("--trim-left" in arguments and not "--trim-right" in arguments) or ("--trim-right" in arguments and not "--trim-left" in arguments):
  trimmed_stats_update(trimmed_stats,sequence,value)
 elif arguments["--operation"]=="trim" and "--trim-right" in arguments and "--trim-left" in arguments :   
  trimmed_stats_update_both(trimmed_stats,sequence,value_left,value_right)
 
#print stats
if arguments["--operation"]=="adaptor-removal" and argument!="--adaptor":
 print("Insert adaptor, please")

if arguments["--operation"]=="trim" and "--trim-left" not in arguments and "--trim-right" not in arguments:
 print("Insert a --trim-left or --trim-right value, please")
  
if arguments["--operation"]=="trim":
 print_summary(arguments,stats,trimmed_stats) 
else:
 print_summary(arguments,stats) 
 
# Close files
input_file.close()
output_file.close()
  


 
 
 
