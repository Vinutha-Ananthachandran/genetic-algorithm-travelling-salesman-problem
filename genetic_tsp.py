import numpy as np
import random

class Location():
    '''
        This class would represent the cities from the input file
        each with their own x,y and z co-ordinate
    '''
    def __init__(self,city):
        '''
            city co-ordinates initialization
        '''
        self.x = city[0]
        self.y = city[1]
        self.z = city[2]
        
    def distance_between_cities(self,city):
        '''
            calculates the distance between the current city object and the city object given in the parameter
            distance formula = sqruareroot of ((x2-x1)^2+(y2-y1)^2+(z2-z1)^2)
        '''
        distance = (((city.x-self.x)**2)+((city.y-self.y)**2)+((city.z-self.z)**2))**0.5
        return distance

class Path():
    '''
        This class contains the path details
    '''
    def __init__(self,path_list):
        self.path_list = path_list
        self.fitness_score = self.fetch_fitness_score()
        
    def fetch_fitness_score(self):
        '''
            traverses the path list and generates the total cost as the fitness score
        '''
        fitness = 0
        #fitness from starting city till the ending city
        for i in range(0,len(self.path_list)-1):
            fitness += self.path_list[i].distance_between_cities(self.path_list[i+1])
        #fitness from ending city back to starting city
        fitness += self.path_list[len(self.path_list)-1].distance_between_cities(self.path_list[0])
        return fitness

class MatingPool():
    '''
        This class contains the functions to generate the mating pool needed for mutation
    '''
    def __init__(self,population):
        self.population = population
        self.rank_list = self.generate_rank_list()
        self.probability_list = self.generate_probability_list()
        self.mating_pool = self.gen_mating_pool()
    
    def gen_mating_pool(self):
        '''
            generates the mating pool using roulette wheel based selection
            and cross breeds the parent paths
        '''
        mating_pool = []
        length = len(self.population)
        i = 0
        while i < length:
            mate = self.roulette_selected_city()
            if mate != -1:
                mating_pool.append(mate)
                i+=1
        return mating_pool

    def generate_rank_list(self):
        '''
            generates the ranklist as a list of tuples where each tuple is of the form
            (path,fitness_score)
        '''
        rank_list = []
        for path in self.population:
            rank_list.append((path,path.fitness_score))
        return rank_list
    
    def generate_probability_list(self):
        '''
            generates the probability list for a given rank list
        '''   
        rank_sum = 0
        for path,rank in self.rank_list:
            rank_sum += rank
        
        self.sort_rank_list()
        
        probability_list = []
        prob_sum = 0
        for path,rank in self.rank_list:
            probability = prob_sum + (rank/rank_sum)
            prob_sum += probability
            probability_list.append((path,probability))
        return probability_list
    
    def sort_rank_list(self):
        '''
            sorts the probabilities in increasing order
        '''
        for i in range(0,len(self.rank_list)):
            for j in range(0,len(self.rank_list)-1-i):
                if self.rank_list[j][1] > self.rank_list[j+1][1]:
                    temp = self.rank_list[j]
                    self.rank_list[j] = self.rank_list[j+1]
                    self.rank_list[j+1] = temp
    
    def roulette_selected_city(self):
        '''
            selects a parent using roulette wheel based selection
        '''            
        rand_probability = np.random.uniform(0,1)
                
        for i in range(0,len(self.probability_list)):
            if self.probability_list[i][1] < rand_probability < self.probability_list[i+1][1]:
                return self.probability_list[i][0]
        return -1

class GenePool():
    '''
        This class contains the functions and attributes related to the gene pool
        which includes cross over and mutation 
    '''
    def __init__(self,mate):
        self.gene_pool = self.generate_gene_pool(mate)
        self.mutated_pool = self.mutate_gene_pool()
        
    def generate_gene_pool(self,mate):
        '''
            fetches 2 parents at a time and cross breeds from parents to generate the gene pool
        '''
        gene = []        
        mlen = len(mate.mating_pool)
        if mlen > 0:
            length = len(mate.mating_pool[0].path_list)
            i = 0
            if mlen <= 1:
                count = mlen
            else:
                count = mlen - 1
            while i < count:
                start = random.randint(0,length-2)
                end = random.randint(start+1,length-1)
                if mlen == 1:
                    parent1,parent2 = mate.mating_pool[i],mate.mating_pool[i]
                else:
                    parent1,parent2 = mate.mating_pool[i],mate.mating_pool[i+1]
                gene.append(self.crossbreed_from_parents(parent1,parent2,start,end))
                i += 1
        return gene
            
    def crossbreed_from_parents(self,parent1,parent2,start,end):
        '''
            fetches elements from start index to end index from parent one 
            and chooses the rest from the parent 2 to create the child
            while following the TSP constraints for each child path
        '''
        #create the child list to be cross bred
        child = [0]*len(parent1.path_list)
        #assign the cities from start index to end index from parent1
        child[start:end+1] = parent1.path_list[start:end+1]
        #find remaining unique cities from parent2 in right of to left order
        index = start
        rem = []
        while len(rem) != (tuple(child).count(0)):
            if index == len(child):
                index = 0
            else:
                if not parent2.path_list[index] in child:
                    rem.append(parent2.path_list[index])
                index += 1
        #populating the child with cities from remaining unique cities
        index = end+1
        while index != start:
            if index == len(child):
                index = 0
            else:
                child[index] = rem.pop(0)
                index += 1
        return Path(child)
    
    def mutate_gene_pool(self):
        '''
            mutates each gene using inversion mutation i.e
            for a given start and stop index, fetch the sub-list and invert(reverse) it
            and place it back into the gene
        '''
        mutated_pool = []
        length = len(self.gene_pool[0].path_list)
        for gene in self.gene_pool:
            start = random.randint(0,length-2)
            end = random.randint(start+1,length-1)
            sub_gene = gene.path_list[start:end+1]
            mutant = gene.path_list[0:start]+sub_gene[::-1]+gene.path_list[end+1::]
            mutated_pool.append(Path(mutant))
        return mutated_pool

class GeneticTravellingSalesman():
    '''
        This class contains the sequence of operations needed to implement travelling sales person problem
    '''
    def __init__(self):
        self.coords = self.gen_input_from_file('input.txt')
        self.cities = []
        if self.coords != None:
            for coord in self.coords:
                self.cities.append(Location(coord))
            self.solution = []
        
    def gen_input_from_file(self,fpath):
        '''
            This function reads the data from the input file;
            Breaks down the data based on the city count in line 1;
            Uses the next n consecutive lines to fetch city coordinates;
        '''
        with open(fpath,mode='r') as file:
            file.seek(0)
            fdata = file.readlines()
        if fdata != []:
            count = int(fdata[0])
            city = []
            for i in range(1,count+1):
                fdata[i] = fdata[i].replace('\n','')
                city.append([int(city) for city in fdata[i].split(' ')])
            return city
    
    def gen_initial_population(self,size):
        '''
            creation of initial population list for parent selection
        '''
        init_pop = []
        for _ in range(0,size):
            rpath = (np.random.permutation(self.cities)).tolist() #random path generated
            if rpath not in init_pop:
                init_pop.append(Path(rpath))
        return init_pop
    
    def travelling_salesman_problem(self):
        '''
            This function employs the culmination of various classes and methods
            in an order to solve TSP
        '''
        best_path = []
        count = len(self.cities)
        if count != 0:
            if count == 1 or count == 2:
                best_path.append(Path(self.cities))
            else:
                if count <= 50:
                    count = count*3 + 200
                elif count > 50 and count <= 100:
                    count += 200
                elif count >=100 and count <=400:
                    count += 100
                else:
                    count = 500
                init_pop = self.gen_initial_population(count)
                for i in range(0,count-1):
                    if i == 0:
                        new_pop = init_pop
                    else:
                        new_pop = gene.mutated_pool
                        del mate
                        del gene
                    mate = MatingPool(new_pop)
                    gene = GenePool(mate)
                    best_path.append(self.fetch_best_path(gene.mutated_pool))
                    del new_pop
            self.gen_output_file(best_path)
        
    def fetch_best_path(self,population):
        '''
            This function fetches the best path in each generation and stores it for 
            optimal final solution
        '''
        min_path = population[0]
        for i in range(1,len(population)):
            if min_path.fitness_score > population[i].fitness_score:
                min_path = population[i]
        return min_path
    
    def gen_output_file(self,best_path):
        '''
            This function generates the output file for valuation
        '''
        self.solution = self.fetch_best_path(best_path)
        with open('output.txt',mode='a') as file:
            for city in self.solution.path_list:
                file.write(f'{city.x} {city.y} {city.z}\n')
            file.write(f'{self.solution.path_list[0].x} {self.solution.path_list[0].y} {self.solution.path_list[0].z}')

tsp = GeneticTravellingSalesman()
tsp.travelling_salesman_problem()