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
#include <random>

using namespace std;

extern __thread std::mt19937 *eng;

class CList;
class Clone;
class MutationHandler;
class OutputWriter;

class ThreadInput{
    /* a datatype for passing in arguments to multiple simulation worker threads.
     should be able to be shared as a single instance for all threads.
     THREAD SAFE for all public methods
     */
private:
    // sim_number is the only member that changes after initialization
    int sim_number;
    // lock for sim_number
    pthread_mutex_t *num_lock;
    string outfolder;
    string infilename;
    string model_type;
    // lock for all writers in the simulation
    pthread_mutex_t *write_lock;
public:
    ThreadInput(pthread_mutex_t *new_lock, pthread_mutex_t *new_write_lock, string new_out, string new_in, string model);
    pthread_mutex_t* getWriteLock(){
        return write_lock;
    }
    /* called by each thread when refreshing the simulation- requests next available simulation number
     @return the ID number of the next simulation. will return even when the next number is greater than the number of simulations to be run.
     */
    int getSimNumberAndAdvance();
    string getOutfolder(){
        return outfolder;
    }
    string getInfile(){
        return infilename;
    }
    string getModel(){
        return model_type;
    }
};

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
    bool has_death_rate;
    double death;
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
    vector<double> empirical_dist;
    vector<double>* getEmpiricalDist(){
        return &empirical_dist;
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
    void addChild(CellType& child_type);
    
    std::vector<CellType *>& getChildren(){
        return children;
    }
    CellType* getParent(){
        return parent;
    }
    void setParent(CellType *parent_type){
        parent = parent_type;
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
    double getDeathRate();

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
    void addDistPoint(double new_birth){
        empirical_dist.push_back(new_birth);
    }
    double getDistByIndex(int index){
        return empirical_dist[index];
    }
    bool hasDist(){
        return empirical_dist.size() > 0;
    }
    int getDistSize(){
        return empirical_dist.size();
    }
    void setDeathRate(double death_rate);
};

class EndListener{
    /* represents a condition under which one simulation trial should be stopped (e.g. after a certain time, after a certain mutation occurs, etc). can have many together in one simulation. do not require refreshing.
     */
public:
    // @return whether the simulation should end
    virtual bool shouldEnd(CList& clone_list) = 0;
    virtual bool readLine(vector<string>& parsed_line) = 0;
};

class MaxCellsListener: public EndListener{
private:
    int max_cells;
public:
    MaxCellsListener();
    bool readLine(vector<string>& parsed_line);
    bool shouldEnd(CList& clone_list);
};

class CompositeListener: public EndListener{
    // container for multiple EndListeners that act on a single simulation
private:
    vector<EndListener*> *listeners;
public:
    CompositeListener();
    ~CompositeListener();
    // @return whether the simulation should end as determined by all of the EndListeners inside
    bool shouldEnd(CList& clone_list);
    bool readLine(vector<string>& parsed_line){return true;}
    void addListener(EndListener& listener);
};

class MaxTimeListener: public EndListener{
private:
    double max_time;
public:
    MaxTimeListener();
    bool readLine(vector<string>& parsed_line);
    bool shouldEnd(CList& clone_list);
};

class OneTypeListener: public EndListener{
public:
    bool readLine(vector<string>& parsed_line){return true;}
    bool shouldEnd(CList& clone_list);
};

class HasTypeListener: public EndListener{
private:
    int type;
    double threshold;
public:
    HasTypeListener();
    bool readLine(vector<string>& parsed_line);
    bool shouldEnd(CList& clone_list);
};

class SimParams{
    /* loads and stores simulation parameters from an input file. see readme for more info on input formatting.
     also initializes key simulation objects (CList and Clones) as appropriate.
     */
private:
    string *model_type;
    int num_simulations;
    int sim_number;
    vector<int> *has_dist;
    vector<vector<int>> *dists;
    MutationHandler *mut_handler;
    CList *clone_list;
    bool handle_line(string& line);
    string *outfolder;
    bool sync_dists;
    bool has_list;
    vector<int> *index_list;
    
    /* handle a line that started with "sim_param".
     @param parsed_line tokenized line with parameter info. already stripped of "sim_param" keyword. first element should be parameter name to be set.
     @return true iff line was correctly formatted and read.
     */
    bool handle_sim_line(vector<string>& parsed_line);
    
    /* all of the following check the parse line, and if possible, makes and inserts the appropriate object into the simulation
     THEY MAY DESTROY OR MODIFY THE PARSED LINE OBJECT
     @return true iff the object was properly created
     */
    bool make_mut_handler();
    bool make_clone(vector<string>& parsed_line);
    bool make_writer(vector<string>& parsed_line);
    bool make_listener(vector<string>& parsed_line);
    bool make_multiclone(vector<string>& parsed_line);
    
    int err_line;
    string err_type;
    string mut_type;
    string sim_name;
    vector<OutputWriter*> *writers;
    CompositeListener *listeners;
    std::vector<string> *mut_params;

public:
    // initial clone list and writer list should be empty, read method will fill it.
    SimParams(CList& clist, vector<OutputWriter*>& writer_list, CompositeListener& listener, string& output, string& sim_type);
    
    /* reads an input file and loads parameters
     MUTATES clone_list, makes new Clones
     @param infile input file to be read
     @return true iff input file was properly formatted and read correctly
     */
    bool read(ifstream& infile);
    void writeErrors(ofstream& errfile);
    int getNumSims(){return num_simulations;}
    string getName(){return sim_name;}
    MutationHandler& get_mut_handler(){return *mut_handler;}
    void refreshSim(ifstream& infile);
    void setSimNumber(int num){
        sim_number = num;
    }
};

#endif /* main_h */
