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
    void beginAction(CList& clone_list){};
};

class CellCountWriter: public DuringOutputWriter{
private:
    int index;
    int sim_number;
    ofstream outfile;
public:
    ~CellCountWriter();
    CellCountWriter(string ofile, int period, int i, int sim);
    void finalAction(CList& clone_list);
    void duringSimAction(CList& clone_list);
    void beginAction(CList& clone_list);
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
    void finalAction(CList& clone_list);
    void duringSimAction(CList& clone_list);
    void beginAction(CList& clone_list);
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
