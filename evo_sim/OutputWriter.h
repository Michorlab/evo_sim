//
//  OutputWriter.hpp
//  evo_sim
//
//  Created by Debra Van Egeren on 2/21/17.
//  Copyright © 2017 Debra Van Egeren. All rights reserved.
//

#ifndef OutputWriter_h
#define OutputWriter_h

#include <stdio.h>
#include <string>
#include <vector>
#include <fstream>

using namespace std;

class CList;

class OutputWriter{
    /* writes results to files. can write/act after every timestep (call to CList::advance()) or only after each simulation instance.
     */
protected:
    string ofile_loc;
    string ofile_name;
    int sim_number;
public:
    virtual void finalAction(CList& clone_list) = 0;
    virtual void duringSimAction(CList& clone_list) = 0;
    virtual void beginAction(CList& clone_list) = 0;
    virtual bool readLine(vector<string>& parsed_line) = 0;
    void setSimNumber(int new_num){
        sim_number = new_num;
    }
    OutputWriter();
};

class FinalOutputWriter: public OutputWriter{
    /* only writes before or after the simulation is run. all writes are thread safe IF AND ONLY IF the write buffer is flushed at the end of each writing method.
     */
public:
    FinalOutputWriter(string ofile);
    virtual void finalAction(CList& clone_list) = 0;
    void duringSimAction(CList& clone_list){};
    virtual void beginAction(CList& clone_list) = 0;
    virtual bool readLine(vector<string>& parsed_line) = 0;
};

class DuringOutputWriter: public OutputWriter{
    /* can update and/or write after every simulation timestep.
     SHOULD NOT write to files shared between multiple trials during the simulation- NOT thread safe. Use a simulation specific file for this purpose.
     */
protected:
    int writing_period;
    int last_written;
    bool shouldWrite(CList& clone_list);
    void resetWriter(){
        last_written = 0;
    };
public:
    DuringOutputWriter(string ofile, int period);
    DuringOutputWriter(string ofile);
    virtual void finalAction(CList& clone_list) = 0;
    virtual void duringSimAction(CList& clone_list) = 0;
    virtual void beginAction(CList& clone_list) = 0;
    virtual bool readLine(vector<string>& parsed_line) = 0;
};

class CountStepWriter: public DuringOutputWriter{
private:
    ofstream outfile;
    int timestep;
    int index;
public:
    CountStepWriter(string ofile);
    ~CountStepWriter();
    bool readLine(vector<string>& parsed_line);
    void finalAction(CList& clone_list);
    void duringSimAction(CList& clone_list);
    void beginAction(CList& clone_list);
    int getTypeIndex(){
        return index;
    }
};

class NumMutationsWriter: public DuringOutputWriter{
private:
    ofstream outfile;
    int index;
public:
    NumMutationsWriter(string ofile);
    ~NumMutationsWriter();
    bool readLine(vector<string>& parsed_line);
    void finalAction(CList& clone_list);
    void duringSimAction(CList& clone_list);
    void beginAction(CList& clone_list);
    int getTypeIndex(){
        return index;
    }
};

class TypeStructureWriter: public FinalOutputWriter{
private:
    ofstream outfile;
public:
    TypeStructureWriter(string ofile);
    ~TypeStructureWriter();
    void finalAction(CList& clone_list);
    void beginAction(CList& clone_list);
    bool readLine(vector<string>& parsed_line){return true;}
};

class CellCountWriter: public DuringOutputWriter{
private:
    int index;
    ofstream outfile;
public:
    ~CellCountWriter();
    CellCountWriter(string ofile, int period, int i, int sim);
    CellCountWriter(string ofile);
    void finalAction(CList& clone_list);
    void duringSimAction(CList& clone_list);
    void beginAction(CList& clone_list);
    bool readLine(vector<string>& parsed_line);
    int getTypeIndex(){
        return index;
    }
};

class FitnessDistWriter: public DuringOutputWriter{
private:
    int index;
    ofstream outfile;
    void write_dist(ofstream& outfile, CList& clone_list);
public:
    ~FitnessDistWriter();
    FitnessDistWriter(string ofile, int period, int i, int sim);
    FitnessDistWriter(string ofile);
    void finalAction(CList& clone_list);
    void duringSimAction(CList& clone_list);
    void beginAction(CList& clone_list);
    bool readLine(vector<string>& parsed_line);
    int getTypeIndex(){
        return index;
    }
};

class MeanFitWriter: public DuringOutputWriter{
private:
    int index;
    ofstream outfile;
public:
    ~MeanFitWriter();
    MeanFitWriter(string ofile, int period, int i, int sim);
    MeanFitWriter(string ofile);
    void finalAction(CList& clone_list);
    void duringSimAction(CList& clone_list);
    void beginAction(CList& clone_list);
    bool readLine(vector<string>& parsed_line);
    int getTypeIndex(){
        return index;
    }
};

class TunnelWriter: public DuringOutputWriter{
private:
    int index;
    bool tunneled;
    ofstream outfile;
public:
    ~TunnelWriter();
    TunnelWriter(string ofile);
    void finalAction(CList& clone_list);
    void duringSimAction(CList& clone_list);
    void beginAction(CList& clone_list);
    bool readLine(vector<string>& parsed_line);
    int getTypeIndex(){
        return index;
    }
};

class AllTypesWriter: public DuringOutputWriter{
private:
    vector<CellCountWriter *> writers;
public:
    ~AllTypesWriter();
    AllTypesWriter(string ofile, int period);
    AllTypesWriter(string ofile);
    void finalAction(CList& clone_list);
    void duringSimAction(CList& clone_list);
    void beginAction(CList& clone_list);
    bool readLine(vector<string>& parsed_line);
};

class IfType2Writer: public FinalOutputWriter{
private:
    ofstream outfile;
public:
    ~IfType2Writer();
    IfType2Writer(string ofile);
    void finalAction(CList& clone_list);
    void beginAction(CList& clone_list){};
    bool readLine(vector<string>& parsed_line){return true;}
};

class IsExtinctWriter: public FinalOutputWriter{
private:
    ofstream outfile;
public:
    ~IsExtinctWriter();
    IsExtinctWriter(string ofile);
    void finalAction(CList& clone_list);
    void beginAction(CList& clone_list){};
    bool readLine(vector<string>& parsed_line){return true;}
};

class EndTimeWriter: public FinalOutputWriter{
private:
    ofstream outfile;
public:
    ~EndTimeWriter();
    EndTimeWriter(string ofile);
    void finalAction(CList& clone_list);
    void beginAction(CList& clone_list){};
    bool readLine(vector<string>& parsed_line){return true;}
};

class EndPopWriter: public FinalOutputWriter{
private:
    ofstream outfile;
public:
    ~EndPopWriter();
    EndPopWriter(string ofile);
    void finalAction(CList& clone_list);
    void beginAction(CList& clone_list){};
    bool readLine(vector<string>& parsed_line){return true;}
};

class NewMutantWriter: public DuringOutputWriter{
private:
    ofstream outfile;
    int index;
    bool has_mutant;
    vector<string> to_write;
public:
    ~NewMutantWriter();
    NewMutantWriter(string ofile);
    void finalAction(CList& clone_list);
    void beginAction(CList& clone_list);
    void duringSimAction(CList& clone_list);
    bool readLine(vector<string>& parsed_line);
};

#endif /* OutputWriter_h */
