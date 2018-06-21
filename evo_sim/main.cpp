#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <chrono>
#include <random>
#include <vector>
#include <unistd.h>
#include <sstream>
#include <string>
#include <algorithm>
#include <pthread.h>


#include "main.h"
#include "Clone.h"
#include "CList.h"
#include "OutputWriter.h"
#include "MutationHandler.h"

// common RNG that is thread safe
__thread std::mt19937 *eng;

void *sim_thread(void *arg){
    int seed1 =  std::chrono::high_resolution_clock::now().time_since_epoch().count();
    eng = new mt19937(seed1);
    ThreadInput *data = (ThreadInput *)arg;
    string outfolder = data->getOutfolder();
    string infilename = data->getInfile();
    string model_type = data->getModel();
    pthread_mutex_t *write_lock = data->getWriteLock();
    
    ofstream errfile;
    errfile.open(outfolder+"input_err.eevo");
    vector<OutputWriter*> writers;
    ifstream infile;
    infile.open(infilename);
    if (!infile.is_open()){
        cout << "bad input file location" << endl;
        cout << infilename << endl;
        pthread_exit(NULL);
    }
    CList* clone_list;
    if (model_type == "moran"){
        clone_list = new MoranPop();
    }
    else if (model_type == "branching"){
        clone_list = new CList();
    }
    else{
        cout << "bad simulation type" << endl;
        pthread_exit(NULL);
    }
    CompositeListener end_conditions;
    SimParams params(*clone_list, writers, end_conditions, outfolder, model_type);
    if (!params.read(infile)){
        params.writeErrors(errfile);
        cout << "bad input file: check error file." << endl;
        infile.close();
        errfile.close();
        pthread_exit(NULL);
    }
    infile.close();
    errfile.close();

    
    int sim_num = data->getSimNumberAndAdvance();
    
    while (sim_num <= params.getNumSims()){
        infile.open(infilename);
        params.refreshSim(infile);
        infile.close();
        params.setSimNumber(sim_num);
        
        pthread_mutex_lock(write_lock);
        for (vector<OutputWriter *>::iterator it = writers.begin(); it != writers.end(); ++it){
            (*it)->setSimNumber(sim_num);
            (*it)->beginAction(*clone_list);
        }
        pthread_mutex_unlock(write_lock);
        while (!clone_list->noTypesLeft() && !clone_list->isExtinct() && !end_conditions.shouldEnd(*clone_list)){
            clone_list->advance();
            for (vector<OutputWriter *>::iterator it = writers.begin(); it != writers.end(); ++it){
                (*it)->duringSimAction(*clone_list);
            }
        }
        pthread_mutex_lock(write_lock);
        for (vector<OutputWriter *>::iterator it = writers.begin(); it != writers.end(); ++it){
            (*it)->finalAction(*clone_list);
        }
        pthread_mutex_unlock(write_lock);
        
        sim_num = data->getSimNumberAndAdvance();
    }
    
    delete clone_list;
    delete eng;
    writers.clear();
    pthread_exit(NULL);
}

int main(int argc, char *argv[]){
    string infilename;
    string outfolder;
    string model_type;
    char tmp;
    int num_cores = 1;
    pthread_mutex_t lock_sim_number;
    pthread_mutex_t lock_writers;
    
    while((tmp=getopt(argc,argv,"i:o:m:n:"))!=-1){
        switch(tmp){
                case 'i':
                infilename = optarg;
                break;
                case 'o':
                outfolder = optarg;
                break;
                case 'm':
                model_type = optarg;
                break;
                case 'n':
                num_cores = stoi(optarg);
                break;
        }
    }
    
    if (infilename=="" || outfolder==""){
        cout << "argument issues" << endl;
        return 1;
    }
    pthread_t threads[num_cores];
    int rc=0;
    pthread_mutex_init(&lock_sim_number, NULL);
    pthread_mutex_init(&lock_writers, NULL);
    ThreadInput thread_data(&lock_sim_number, &lock_writers, outfolder, infilename, model_type);
    
    for (int i=0; i<num_cores; i++){
        
        if (rc != pthread_create(&threads[i], NULL, sim_thread, &thread_data)){
            cout << "thread creation failure" << endl;
            return 1;
        }
    }
    for (int i = 0; i < num_cores; ++i) {
        pthread_join(threads[i], NULL);
    }
    return 0;
}

//=============CLASS METHODS==================

ThreadInput::ThreadInput(pthread_mutex_t *new_lock, pthread_mutex_t *new_write_lock, string new_out, string new_in, string model){
    sim_number = 1;
    num_lock = new_lock;
    write_lock = new_write_lock;
    outfolder = new_out;
    infilename = new_in;
    model_type = model;
}

int ThreadInput::getSimNumberAndAdvance(){
    pthread_mutex_lock(num_lock);
    int sim_num = sim_number;
    sim_number++;
    pthread_mutex_unlock(num_lock);
    return sim_num;
}

CellType::CellType(int i, CellType *parent_type){
    index = i;
    total_birth_rate = 0;
    parent = parent_type;
    children = std::vector<CellType *>();
    if (parent_type){
        empirical_dist = *parent_type->getEmpiricalDist();
    }
    num_cells = 0;
    root_node = NULL;
    end_node = NULL;
    prev_node = NULL;
    next_node = NULL;
    has_death_rate = false;
    death = 0.0;
}

void CellType::setDeathRate(double death_rate){
    has_death_rate = true;
    death = death_rate;
    clone_list->death_var = true;
}

double CellType::getDeathRate(){
    if (has_death_rate){
        return death;
    }
    else{
        return clone_list->getDeathRate();
    }
}

void CellType::subtractOneCell(double b){
    num_cells --;
    total_birth_rate-=b;
    /*
    if (isExtinct()){
        unlinkType();
    }
     */
    clone_list->removeCell(b);
}

void CellType::addChild(CellType &child_type){
    bool has_child = false;
    for (vector<CellType *>::iterator it = children.begin(); it != children.end(); ++it){
        if ((*it)->getIndex() == child_type.getIndex()){
            has_child = true;
        }
    }
    if (!has_child){
        children.push_back(&child_type);
    }
}

void CellType::addCells(int num, double b){
    num_cells += num;
    total_birth_rate += b*num;
    
    clone_list->addCells(num, b);
}

CellType::~CellType(){
    Clone *to_delete = root_node;
    Clone *next;
    while (to_delete){
        next = &to_delete->getNextWithinType();
        delete to_delete;
        to_delete = next;
    }
}

void CellType::unlinkType(){
    if (next_node){
        next_node->setPrev(*prev_node);
    }
    if (prev_node){
        prev_node->setNext(*next_node);
    }
}

MutationHandler& CellType::getMutHandler(){
    return clone_list->getMutHandler();
}

void CellType::insertClone(Clone &new_clone){
    addCells(new_clone.getCellCount(), new_clone.getBirthRate());
    if (!root_node){
        root_node = &new_clone;
        end_node = &new_clone;
    }
    else{
        end_node->setNext(&new_clone);
        new_clone.setPrev(end_node);
        end_node = &new_clone;
    }
}

//----------EndListeners----------------
MaxTimeListener::MaxTimeListener(){
    max_time = 0;
}

MaxCellsListener::MaxCellsListener(){
    max_cells = 0;
}

bool MaxCellsListener::readLine(vector<string>& parsed_line){
    try {
        max_cells =stoi(parsed_line[0]);
    }
    catch (...){
        return false;
    }
    return true;
}

bool MaxTimeListener::readLine(vector<string>& parsed_line){
    try {
        max_time = stod(parsed_line[0]);
    }
    catch (...){
        return false;
    }
    return true;
}

bool MaxCellsListener::shouldEnd(CList& clone_list){
    return clone_list.getNumCells() >= max_cells;
}

bool MaxTimeListener::shouldEnd(CList& clone_list){
    return clone_list.getCurrTime() >= max_time;
}

bool HasTypeListener::readLine(vector<string>& parsed_line){
    try {
        type =stoi(parsed_line[0]);
        threshold =stod(parsed_line[1]);
    }
    catch (...){
        return false;
    }
    return true;
}

HasTypeListener::HasTypeListener(){
    type = 0;
    threshold = 0;
}

bool HasTypeListener::shouldEnd(CList& clone_list){
    CellType *new_type = clone_list.getTypeByIndex(type);
    return new_type && new_type->getNumCells() >= threshold;
}

bool OneTypeListener::shouldEnd(CList& clone_list){
    return clone_list.isOneType();
}

CompositeListener::CompositeListener(){
    listeners = new vector<EndListener*>();
}

bool CompositeListener::shouldEnd(CList& clone_list){
    bool flag = false;
    for (vector<EndListener *>::iterator it = listeners->begin(); it != listeners->end(); ++it){
        flag = flag || (*it)->shouldEnd(clone_list);
    }
    return flag;
}

void CompositeListener::addListener(EndListener& listener){
    listeners->push_back(&listener);
}

CompositeListener::~CompositeListener(){
    listeners->clear();
}

//--------------SimParams-----------

SimParams::SimParams(CList& clist, vector<OutputWriter*>& writer_list, CompositeListener& listener, string& output, string& sim_type){
    num_simulations = 0;
    mut_handler = NULL;
    err_type = "";
    err_line = 0;
    mut_type = "";
    mut_params = NULL;
    sim_name = "";
    clone_list = &clist;
    sim_number = 1;
    writers = &writer_list;
    outfolder = &output;
    listeners = &listener;
    has_dist = new vector<int>();
    dists = new vector<vector<int>>();
    sync_dists = false;
    has_list = false;
    index_list = new vector<int>();
    model_type = &sim_type;
}

void SimParams::refreshSim(ifstream& infile){
    clone_list->refreshSim();
    string line;
    
    int num_drawn = -1;
    string tok;
    std::vector<string> parsed_line;
    
    while (getline(infile, line)){
        std::stringstream ss;
        ss.str(line);
        while (getline(ss, tok, ' ')){
            parsed_line.push_back(tok);
        }
        if (parsed_line[0] == "clone"){
            parsed_line.erase(parsed_line.begin());
            int index = stoi(parsed_line[1]);
            auto found = std::find(has_dist->begin(), has_dist->end(), index);
            if (found != has_dist->end()){
                auto new_index = std::distance(has_dist->begin(), found);
                int num_cells = 0;
                if (has_list){
                    if (sim_number >= index_list->size()){
                        cout << "index list too short for specified number of simulations" << endl;
                    }
                    int list_index = index_list->at(sim_number-1);
                    if (list_index < 0 || list_index > dists->at(new_index).size()){
                        cout << "list contains bad indices" << endl;
                    }
                    num_cells = dists->at(new_index).at(list_index);
                }
                else if (sync_dists && num_drawn >= 0){
                    if (num_drawn > dists->at(new_index).size()){
                        cout << "synched cell count dists don't have same size" << endl;
                    }
                    num_cells = dists->at(new_index).at(num_drawn);
                }
                else{
                    uniform_real_distribution<double> runif;
                    int ran = floor(runif(*eng) * dists->at(new_index).size());
                    num_cells = dists->at(new_index).at(ran);
                    num_drawn = ran;
                    
                }
                parsed_line[2] = to_string(num_cells);
            }
            make_clone(parsed_line);
        }
        else if(parsed_line[0] == "multiclone"){
            parsed_line.erase(parsed_line.begin());
            make_multiclone(parsed_line);
        }
        parsed_line.clear();
    }

}

bool SimParams::read(ifstream& infile){
    string line;
    
    int line_num = 1;
    while (getline(infile, line)){
        if (!handle_line(line)){
            err_line = line_num;
            return false;
        }
        line_num++;
    }
    if (!make_mut_handler()){
        return false;
    }
    else if (clone_list->isExtinct()){
        err_type = "no clones added successfully";
        return false;
    }
    clone_list->setMutHandler(*mut_handler);
    return true;
}

bool SimParams::handle_line(string& line){
    std::stringstream ss;
    ss.str(line);
    string tok;
    std::vector<string> parsed_line;
    while (getline(ss, tok, ' ')){
        parsed_line.push_back(tok);
    }
    if (parsed_line[0][0] == '#' || parsed_line.size()==0){
        return true;
    }
    else if (parsed_line[0] == "sim_params"){
        parsed_line.erase(parsed_line.begin());
        handle_sim_line(parsed_line);
    }
    else if (parsed_line[0] == "pop_params"){
        parsed_line.erase(parsed_line.begin());
        if (!clone_list->handle_line(parsed_line)){
            err_type = "bad pop_params line";
            return false;
        }
    }
    else if (parsed_line[0] == "init_dist"){
        parsed_line.erase(parsed_line.begin());
        if (parsed_line[0] == "sync"){
            sync_dists = true;
            return true;
        }
        else if (parsed_line.size() != 2){
            return false;
        }
        else if (parsed_line[0] == "list"){
            if (has_list){
                err_type = "two init_dist list lines";
                return false;
            }
            has_list = true;
            ifstream infile;
            infile.open(parsed_line[1]);
            if (!infile.is_open()){
                err_type = "bad init_dist list file";
                return false;
            }
            string line;
            while (getline(infile, line)){
                std::stringstream ss;
                string tok;
                ss.str(line);
                try{
                    while (getline(ss, tok)){
                        index_list->push_back(stod(tok));
                    }
                }
                catch(...){
                    return false;
                }
            }
            return true;

        }
        has_dist->push_back(stoi(parsed_line[0]));
        dists->push_back(*new vector<int>());
        ifstream infile;
        infile.open(parsed_line[1]);
        if (!infile.is_open()){
            err_type = "bad init_dist file";
            return false;
        }
        string line;
        while (getline(infile, line)){
            std::stringstream ss;
            string tok;
            ss.str(line);
            try{
                while (getline(ss, tok)){
                    dists->back().push_back(stod(tok));
                }
            }
            catch(...){
                return false;
            }
        }

    }
    else if (parsed_line[0] == "clone"){
        parsed_line.erase(parsed_line.begin());
        if (!make_clone(parsed_line)){
            return false;
        }
    }
    else if (parsed_line[0] == "multiclone"){
        /* used when creating several clone types which are equivalent except in type id
         syntax: multiclone [num_types] clone [clone_type] [clone params]*
         the type id provided in the clone params IS THE LOWEST CLONE ID OF THE SET OF CLONES TO BE MADE
         */
        parsed_line.erase(parsed_line.begin());
        if (!make_multiclone(parsed_line)){
            return false;
        }
    }
    else if (parsed_line[0] == "writer"){
        parsed_line.erase(parsed_line.begin());
        if (!make_writer(parsed_line)){
            return false;
        }
    }
    else if (parsed_line[0] == "listener"){
        parsed_line.erase(parsed_line.begin());
        if (!make_listener(parsed_line)){
            return false;
        }
    }
    else{
        err_type = "bad first keyword " + parsed_line[0];
        return false;
    }
    return true;
}

bool SimParams::make_multiclone(vector<string> &parsed_line){
    int num_types = stoi(parsed_line[0]);
    parsed_line.erase(parsed_line.begin());
    parsed_line.erase(parsed_line.begin());
    std::vector<string> parsed_copy = parsed_line;
    for (int i=0; i<num_types; i++){
        parsed_line = parsed_copy;
        parsed_line[1] = to_string(i + stoi(parsed_line[1]));
        if (!make_clone(parsed_line)){
            return false;
        }
    }
    return true;
}

bool SimParams::make_listener(vector<string> &parsed_line){
    if (parsed_line.size() < 1){
        err_type = "bad listener parameters";
        return false;
    }
    string type = parsed_line[0];
    parsed_line.erase(parsed_line.begin());
    EndListener *new_listener;
    if (type == "HasType"){
        new_listener = new HasTypeListener();
    }
    else if (type == "OneType"){
        new_listener = new OneTypeListener();
        
    }
    else if (type == "MaxTime"){
        new_listener = new MaxTimeListener();
    }
    else if (type == "MaxCells"){
        new_listener = new MaxCellsListener();
    }
    else{
        err_type = "bad listener type";
        return false;
    }
    if (!new_listener->readLine(parsed_line)){
        err_type = "bad listener params";
        return false;
    }
    listeners->addListener(*new_listener);
    return true;
}

bool SimParams::make_writer(vector<string> &parsed_line){
    if (parsed_line.size() < 1){
        err_type = "bad writer parameters";
        return false;
    }
    OutputWriter *new_writer;
    string type = parsed_line[0];
    parsed_line.erase(parsed_line.begin());
    if (type == "IfType2"){
        new_writer = new IfType2Writer(*outfolder);
    }
    else if (type == "IsExtinct"){
        new_writer = new IsExtinctWriter(*outfolder);
    }
    else if (type == "IfType"){
        new_writer = new IfTypeWriter(*outfolder);
    }
    else if (type == "TypeStructure"){
        new_writer = new TypeStructureWriter(*outfolder);
    }
    else if (type == "AllTypes"){
        new_writer = new AllTypesWriter(*outfolder);
    }
    else if (type == "CellCount"){
        new_writer = new CellCountWriter(*outfolder);
    }
    else if (type == "FitnessDist"){
        new_writer = new FitnessDistWriter(*outfolder);
    }
    else if (type == "MotherDaughter"){
        new_writer = new MotherDaughterWriter(*outfolder);
    }
    else if (type == "MeanFit"){
        new_writer = new MeanFitWriter(*outfolder);
    }
    else if (type == "NewMutant"){
        new_writer = new NewMutantWriter(*outfolder);
    }
    else if (type == "EndTime"){
        new_writer = new EndTimeWriter(*outfolder);
    }
    else if (type == "EndPop"){
        new_writer = new EndPopWriter(*outfolder);
    }
    else if (type == "EndPopTypes"){
        new_writer = new EndPopTypesWriter(*outfolder);
    }
    else if (type == "CountStep"){
        new_writer = new CountStepWriter(*outfolder);
    }
    else if (type == "Tunnel"){
        new_writer = new TunnelWriter(*outfolder);
    }
    else if (type == "NumMutations"){
        new_writer = new NumMutationsWriter(*outfolder);
    }
    else{
        err_type = "bad writer type";
        return false;
    }
    if (!new_writer->readLine(parsed_line)){
        err_type = "bad writer params";
        return false;
    }
    writers->push_back(new_writer);
    return true;
}

bool SimParams::make_clone(vector<string> &parsed_line){
    if (parsed_line.size() < 3){
        err_type = "bad params for Clone";
        return false;
    }
    string type = parsed_line[0];
    parsed_line.erase(parsed_line.begin());
    int type_id =stoi(parsed_line[0]);
    int num_cells =stoi(parsed_line[1]);
    
    if (clone_list->getTypeByIndex(type_id)){
        err_type = "type space conflict";
        return false;
    }
    
    CellType *new_type = new CellType(type_id, NULL);
    
    clone_list->addRootType(*new_type);
    clone_list->insertCellType(*new_type);
    parsed_line.erase(parsed_line.begin());
    if (type == "Simple"){
        if (parsed_line.size() < 3){
            err_type = "bad params for SimpleClone";
            return false;
        }
        Clone *new_clone = new SimpleClone(*new_type);
        if (!new_clone->readLine(parsed_line)){
            return false;
        }
        new_type->insertClone(*new_clone);
    }
    else if (type == "TypeSpecific"){
        if (parsed_line.size() < 4){
            err_type = "bad params for TypeSpecificClone";
            return false;
        }
        TypeSpecificClone *new_clone;
        for (int i=0; i<num_cells; i++){
            if (*model_type == "moran"){
                new_clone = new TypeSpecificClone(*new_type, true);
            }
            else{
                new_clone = new TypeSpecificClone(*new_type, false);
            }
            if (!new_clone->readLine(parsed_line)){
                return false;
            }
            new_type->insertClone(*new_clone);
        }
    }
    else if (type == "Heritable"){
        if (parsed_line.size() < 4){
            err_type = "bad params for HeritableClone";
            return false;
        }
        Clone *new_clone;
        for (int i=0; i<num_cells; i++){
            
            if (*model_type == "moran"){
                new_clone = new HeritableClone(*new_type, true);
            }
            else{
                new_clone = new HeritableClone(*new_type, false);
            }
            if (!new_clone->readLine(parsed_line)){
                err_type = "bad clone";
                return false;
            }
            new_type->insertClone(*new_clone);
        }
    }
    else if (type == "HerReset"){
        if (parsed_line.size() < 4){
            err_type = "bad params for HerResetClone";
            return false;
        }
        Clone *new_clone;
        for (int i=0; i<num_cells; i++){
            
            if (*model_type == "moran"){
                new_clone = new HerResetClone(*new_type, true);
            }
            else{
                new_clone = new HerResetClone(*new_type, false);
            }
            if (!new_clone->readLine(parsed_line)){
                err_type = "bad clone";
                return false;
            }
            new_type->insertClone(*new_clone);
        }
    }
    else if (type == "HerResetExp"){
        if (parsed_line.size() < 5){
            err_type = "bad params for HerResetExpClone";
            return false;
        }
        Clone *new_clone;
        for (int i=0; i<num_cells; i++){
            
            if (*model_type == "moran"){
                new_clone = new HerResetExpClone(*new_type, true);
            }
            else{
                new_clone = new HerResetExpClone(*new_type, false);
            }
            if (!new_clone->readLine(parsed_line)){
                err_type = "bad clone";
                return false;
            }
            new_type->insertClone(*new_clone);
        }
    }
    else if (type == "HerPoisson"){
        if (parsed_line.size() < 4){
            err_type = "bad params for HerPoissonClone";
            return false;
        }
        Clone *new_clone;
        for (int i=0; i<num_cells; i++){
            
            if (*model_type == "moran"){
                new_clone = new HerPoissonClone(*new_type, true);
            }
            else{
                new_clone = new HerPoissonClone(*new_type, false);
            }
            if (!new_clone->readLine(parsed_line)){
                err_type = "bad clone";
                return false;
            }
            new_type->insertClone(*new_clone);
        }
    }
    else if (type == "HerResetEmpiric"){
        if (parsed_line.size() < 5){
            err_type = "bad params for HerResetEmpiricClone";
            return false;
        }
        Clone *new_clone;
        for (int i=0; i<num_cells; i++){
            
            if (*model_type == "moran"){
                new_clone = new HerResetEmpiricClone(*new_type, true);
            }
            else{
                new_clone = new HerResetEmpiricClone(*new_type, false);
            }
            if (!new_clone->readLine(parsed_line)){
                err_type = "bad clone";
                return false;
            }
            new_type->insertClone(*new_clone);
        }
    }
    else if (type == "HerEmpiric"){
        if (parsed_line.size() < 5){
            err_type = "bad params for HerEmpiricClone";
            return false;
        }
        Clone *new_clone;
        for (int i=0; i<num_cells; i++){
            if (*model_type == "moran"){
                new_clone = new HerEmpiricClone(*new_type, true);
            }
            else{
                new_clone = new HerEmpiricClone(*new_type, false);
            }
            if (!new_clone->readLine(parsed_line)){
                err_type = "bad clone";
                return false;
            }
            new_type->insertClone(*new_clone);
        }
    }
    else if (type == "TypeEmpiric"){
        if (parsed_line.size() < 5){
            err_type = "bad params for TypeEmpiricClone";
            return false;
        }
        Clone *new_clone;
        for (int i=0; i<num_cells; i++){
            if (*model_type == "moran"){
                new_clone = new TypeEmpiricClone(*new_type, true);
            }
            else{
                new_clone = new TypeEmpiricClone(*new_type, false);
            }
            
            if (!new_clone->readLine(parsed_line)){
                err_type = "bad clone";
                return false;
            }
            new_type->insertClone(*new_clone);
        }
    }
    else{
        err_type = "bad clone type information";
        return false;
    }
    return true;
}

bool SimParams::handle_sim_line(vector<string>& parsed_line){
    
    if (parsed_line.size() == 1 && parsed_line[0].at(0) == '#'){
        return true;
    }
    else if (parsed_line.size() <= 1){
        return false;
    }
    if (parsed_line[0] == "num_simulations"){
        num_simulations =stoi(parsed_line[1]);
    }
    else if (parsed_line[0] == "mut_handler_type"){
        mut_type = parsed_line[1];
    }
    else if (parsed_line[0] == "mut_handler_params"){
        parsed_line.erase(parsed_line.begin());
        mut_params = new std::vector<string>(parsed_line);
    }
    else if (parsed_line[0] == "sim_id"){
        sim_name = parsed_line[1];
    }
    return true;
}

bool SimParams::make_mut_handler(){
    if (mut_type == "ThreeTypes"){
        mut_handler = new ThreeTypesMutation();
    }
    else if (mut_type == "ThreeTypesMult"){
        mut_handler = new ThreeTypesMultMutation();
    }
    else if (mut_type == "FixedSites"){
        mut_handler = new FixedSitesMutation();
    }
    else if (mut_type == "Neutral"){
        mut_handler = new NeutralMutation();
    }
    else if (mut_type == "None"){
        mut_handler = new NoMutation();
    }
    else if (mut_type == "ThreeTypesFlex"){
        mut_handler = new ThreeTypesFlexMutation();
    }
    else if (mut_type == "ManyTypesFlex"){
        mut_handler = new ManyTypesFlexMutation();
    }
    else{
        err_type = "bad mut type";
        return false;
    }
    if (!mut_handler->read(*mut_params)){
        err_type = "bad mut params";
        return false;
    }
    return true;
}

void SimParams::writeErrors(ofstream& errfile){
    errfile << "SIM INPUT ERRORS" << endl;
    errfile << "error type: " << err_type << endl;
    errfile << "error line: " << err_line << endl;
}
