#ifndef main_h
#define main_h

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <vector>

using namespace std;

class CList;
class Clone;
class MutationHandler;
class OutputWriter;

class CellType{
    /* represents a functional subset of cells in the population (e.g. cells with a specific mutation, phenotype, etc)
     distinct from fitness- cells with different birth rates can have the same type
     tracks the type history/phylogeny and number of cells by type
     */
    friend class CList;
    friend class Clone;
private:
    CellType *parent;
    std::vector<CellType *> children;
    CellType *prev_node;
    CellType *next_node;
    Clone *root_node;
    Clone *end_node;
    int index;
    int num_cells;
    double total_birth_rate;
    void unlinkType();
    void setNext(CellType& next){
        next_node = &next;
    }
    void setPrev(CellType& prev){
        prev_node = &prev;
    }
    CList *clone_list;
    void setRoot(Clone& new_root){
        root_node = &new_root;
    }
    void setEnd(Clone& new_end){
        end_node = &new_end;
    }
    // called every time a cell of this type is born
    void addCells(int num, double b);
    // called every time a cell of this type dies
    void subtractOneCell(double b);
    void setCloneList(CList& clist){
        clone_list = &clist;
    }
public:
    /* @param i cell type id. should be unique in the clone list typespace.
     @param parent_type cell type that formed this type, via mutation. if one of the original types in simulation, then NULL.
     @param num initial number of cells of this type.
     */
    CellType(int i, CellType *parent_type);
    
    ~CellType();
    
    /* called when a new type is formed after mutation from this parent type
     @param child_type child to be added
     */
    void addChild(CellType& child_type){
        children.push_back(&child_type);
    }
    std::vector<CellType *>& getChildren(){
        return children;
    }
    CellType* getParent(){
        return parent;
    }
    bool isExtinct(){
        return num_cells == 0;
    }
    int getIndex(){
        return index;
    }
    
    int getNumCells(){
        return num_cells;
    }
    double getBirthRate(){
        return total_birth_rate;
    }

    CellType* getNext(){
        return next_node;
    }
    Clone* getRoot(){
        return root_node;
    }
    Clone* getEnd(){
        return end_node;
    }
    MutationHandler& getMutHandler();
    CList& getPopulation(){
        return *clone_list;
    }
    void insertClone(Clone& new_clone);
};

class MutationHandler {
    /* chooses new parameters (mutation prob, fitness, cell type) for a new cell type
     use: first call generateMutant, then get new parameters
     ABSTRACT CLASS
     */
protected:
    double birth_rate;
    double mut_prob;
    CellType *new_type;
    
    /* gets the desired mutant type, either the existing type or creates a new type. 
     @param index desired index of new type. SHOULD NOT be outside the range 0<=index<max_types of clone_list
     @param clone_list the population in which the new type will ultimately reside
     @param curr_type the parent of the new cell
     @return a cell type appropriate for the typespace in clone_list with the desired index. MUST BE DELETED LATER.
     */
    CellType* getNewTypeByIndex(int index, CellType& curr_type);
public:
    //@return birth rate of new type, after generateMutant
    double getNewBirthRate() {return birth_rate;}
    
    //@return mutation probability of new type, after generateMutant
    double getNewMutProb() {return mut_prob;}
    
    //@return CellType of new type, after generateMutant. should not conflict with current population typespace in clone_list.
    CellType& getNewType() {return *new_type;}
    
    /*
     loads parameters and calculates mutant daughter mutation rate, type, and birth rate
     @param type current type of mother cell
     @param next_type unused type id that can be taken by a new mutant if needed
     @param b current birth rate of mother
     @param mut current mutation rate of mother
     */
    virtual void generateMutant(CellType& type, double b, double mut) = 0;
    virtual bool read(std::vector<string>& params) = 0;
};

class ThreeTypesMutation: public MutationHandler {
    /* only forward mutation, three types {0,1,2} in typespace, constant mutation rates between types
     additive fitness changes between types
     conceptually similar to 2-hit TS model
     type 2 absorbing state
     */
private:
    double mu2;
    double fit1;
    double fit2;
public:
    ThreeTypesMutation(){};
    ThreeTypesMutation(double m2, double f1, double f2);
    void generateMutant(CellType& type, double b, double mut);
    bool read(std::vector<string>& params);
};

class SimParams{
    /* loads and stores simulation parameters from an input file. see readme for more info on input formatting.
     also initializes key simulation objects (CList and Clones) as appropriate.
     */
private:
    int num_simulations;
    double max_time;
    int max_cells;
    MutationHandler *mut_handler;
    CList *clone_list;
    bool handle_line(string& line);
    string *outfolder;
    
    /* handle a line that started with "sim_param".
     @param parsed_line tokenized line with parameter info. already stripped of "sim_param" keyword. first element should be parameter name to be set.
     @return true iff line was correctly formatted and read.
     */
    bool handle_sim_line(vector<string>& parsed_line);
    bool make_mut_handler();
    bool make_clone(vector<string>& parsed_line);
    bool make_writer(vector<string>& parsed_line);
    int err_line;
    string err_type;
    string mut_type;
    string sim_name;
    vector<OutputWriter*> *writers;
    std::vector<string> *mut_params;

public:
    // initial clone list and writer list should be empty, read method will fill it.
    SimParams(CList& clist, vector<OutputWriter*>& writer_list, string& output);
    
    /* reads an input file and loads parameters
     MUTATES clone_list, makes new Clones
     @param infile input file to be read
     @return true iff input file was properly formatted and read correctly
     */
    bool read(ifstream& infile);
    void writeErrors(ofstream& errfile);
    int getNumSims(){return num_simulations;}
    double getMaxTime(){return max_time;}
    int getMaxCells(){return max_cells;}
    string getName(){return sim_name;}
    MutationHandler& get_mut_handler(){return *mut_handler;}
    void refreshSim(ifstream& infile);
};

#endif /* main_h */
