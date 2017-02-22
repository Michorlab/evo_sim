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
public:
    virtual void finalAction(CList& clone_list) = 0;
    virtual void duringSimAction(CList& clone_list) = 0;
    virtual void beginAction(CList& clone_list) = 0;
};

class FinalOutputWriter: public OutputWriter{
protected:
    int sim_number;
public:
    FinalOutputWriter(string ofile);
    virtual void finalAction(CList& clone_list) = 0;
    void duringSimAction(CList& clone_list){};
    virtual void beginAction(CList& clone_list) = 0;
};

class DuringOutputWriter: public OutputWriter{
protected:
    int writing_period;
    int last_written;
    bool shouldWrite(CList& clone_list);
public:
    DuringOutputWriter(string ofile, int period);
    DuringOutputWriter(string ofile);
    virtual void finalAction(CList& clone_list) = 0;
    virtual void duringSimAction(CList& clone_list) = 0;
    virtual void beginAction(CList& clone_list) = 0;
};

class TypeStructureWriter: public FinalOutputWriter{
private:
    ofstream outfile;
public:
    TypeStructureWriter(string ofile);
    ~TypeStructureWriter();
    void finalAction(CList& clone_list);
    void beginAction(CList& clone_list);
};

class CellCountWriter: public DuringOutputWriter{
private:
    int index;
    int sim_number;
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
    int sim_number;
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

class AllTypesWriter: public DuringOutputWriter{
private:
    int sim_number;
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
};

class IsExtinctWriter: public FinalOutputWriter{
private:
    ofstream outfile;
public:
    ~IsExtinctWriter();
    IsExtinctWriter(string ofile);
    void finalAction(CList& clone_list);
    void beginAction(CList& clone_list){};
};

class EndTimeWriter: public FinalOutputWriter{
private:
    ofstream outfile;
public:
    ~EndTimeWriter();
    EndTimeWriter(string ofile);
    void finalAction(CList& clone_list);
    void beginAction(CList& clone_list){};
};


#endif /* OutputWriter_h */