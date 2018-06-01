// based off https://github.com/mokarrom/mpi-vina

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>

#include "mpi.h"
// #include "CLinkedList.h"

#define MASTER                  0

#define COMPUTE_TAG             11
#define TERMINATE_TAG           22
#define WORK_REQ_TAG            33

#define MAX_LIGAND_NAME_LENGTH  25
#define LIGAND_FILE_NAME        "ligandlist"
#define TMPDIR                  "/tmp/"
// TODO: 
// pipe output to database

/*
 * It will iterate through all the lines in file and
 * put them in given vector
 */
bool getFileContent(std::string fileName, std::vector<std::string> & vecOfStrs)
{
        // Open the File
        std::ifstream in(fileName.c_str());

        // Check if object is valid
        if(!in)
        {
                printf("Cannot open the file %s.\n", fileName.c_str());
                return false;
        }

        std::string str;
        // Read the next line from File untill it reaches the end.
        while (std::getline(in, str))
        {
                // Line contains string of length > 0 then save it in vector
                if(str.size() > 0)
                        vecOfStrs.push_back(str);
        }
        //Close The File
        in.close();
        return true;
}

void sminaMpiManager (int numProcs);
void sminaMpiWorker (int workerId, int coresPerTask);

//MPI_Datatype MPI_LIGAND;
// LigandList lgndsList;
std::vector<std::string> ligandNames;

int main(int argc, char *argv[])
{
    if(argc!=3)
    {
        printf("Usage: %s <ligandlist> <CoresPerTask>\n\tligandlist: list of ligand names, one per line\n\tCoresPerTask: cpu cores used for each smina task\n", "sminaMPI");
        return 1;
    }

    int numProcs, rank;
    double startTime = 0.0, endTime = 0.0;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    //MPI_Type_contiguous(MAX_LIGAND_NAME_LENGTH, MPI_CHAR, &MPI_LIGAND);
    //MPI_Type_commit(&MPI_LIGAND);

    int totalLigands;
    int coresPerTask = 1;
    coresPerTask = atoi(argv[2]);
    // Only master processor will read the ligandlist file and will make the work pool.
    if(rank == MASTER)
    {
        printf("Master processor : Reading ligandlist file and creating work pool.\n");
        startTime = MPI_Wtime(); // start timer.

        //ligandNames.reserve(1000000);
	bool getContentOK = getFileContent(argv[1], ligandNames);

        if (! getContentOK)
        {
            printf("Couldn't open file %s for reading.\n", argv[1]);
            MPI_Abort(MPI_COMM_WORLD, 911); //Terminates MPI execution environment with error code 911.
            return 1;
        }

        //Get total number of ligands.
        totalLigands = ligandNames.size();
        printf("Total ligands loaded from ligandlist(%s): %u in %.2lf seconds. front=%s back=%s\n", argv[1], totalLigands, MPI_Wtime() - startTime, ligandNames.front().c_str(), ligandNames.back().c_str()); 
        fflush(stdout);
    }

    if (rank == MASTER)
    {
        sminaMpiManager(numProcs);   // Master processor will play the role of smina manager.
    }
    else
    {
        sminaMpiWorker(rank, coresPerTask);    // All other processors will play the role of smina worker.
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == MASTER)
    {
        endTime = MPI_Wtime(); // end timer.
        printf("\n\n..........................................\n"); 
        printf("   Number of workers       = %d \n", numProcs - 1); 
        printf("   Number of Lignds        = %u \n", totalLigands); 
        printf("   Total time required     = %.2lf seconds.\n", endTime - startTime); 
        printf("..........................................\n\n");
		fflush(stdout);

        ligandNames.clear();
        std::vector<std::string>().swap(ligandNames); //free up mem
    }

    //MPI_Type_free (&MPI_LIGAND);
    MPI_Finalize ( );
    return 0;
}

void sminaMpiManager(int numProcs)
{
    int i = 0;
    MPI_Status mStatus;
    time_t tStart  = time(NULL);

    // for progress
    const int percentPrint = 5;
    int step = ligandNames.size() / (100/percentPrint);
    int nextPrint = step;

    //std::string currLigand = ligandNames[i];    //Begin with first ligand.
    //MPI_Barrier(MPI_COMM_WORLD);
    printf("Master has started.\n");
	fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    while (i < ligandNames.size())
    {
        //Wait for receiving work item request from any worker.
        MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, WORK_REQ_TAG, MPI_COMM_WORLD, &mStatus);
        //Assign an work item to the requested worker.
        MPI_Send(ligandNames[i].c_str(), ligandNames[i].length(), MPI_CHAR, mStatus.MPI_SOURCE, COMPUTE_TAG, MPI_COMM_WORLD);
        //printf("INFO: send %d char: %s.\n", ligandNames[i].length(), ligandNames[i].c_str());

        // progress
        if (i+1 >= nextPrint)
        {
            int percent = (100 * (i+1)) / ligandNames.size();
            printf("INFO: send jobs %d / %d (%d%%) in %d seconds.\n", i+1, ligandNames.size(), percent, time(NULL) - tStart);
            fflush(stdout);
            nextPrint += step;
        }

        i++;    //Go ahead of the list.
    }

    //Computation has done. Send termination tag to all the slaves.
    for (unsigned int j = 1; j < numProcs; j++)
    {
        MPI_Recv(NULL, 0, MPI_INT, MPI_ANY_SOURCE, WORK_REQ_TAG, MPI_COMM_WORLD, &mStatus);
        MPI_Send(NULL, 0, MPI_INT, mStatus.MPI_SOURCE, TERMINATE_TAG, MPI_COMM_WORLD);
    }
    printf("INFO: all TERMINATE_TAG sent. sminaMpiManager exiting\n");
    //MPI_Barrier(MPI_COMM_WORLD);
    return;
}

void sminaMpiWorker(int workerId, int coresPerTask)
{
    MPI_Status wStatus;
    char thisLigandName[MAX_LIGAND_NAME_LENGTH];
    int nameLen;
    int numJobs = 0;
    double startTime = 0.0, endTime = 0.0;
    startTime = MPI_Wtime(); // start timer.

    // Get the number of processes and the name of the processor
    int world_size;
    int processor_name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Get_processor_name(processor_name, &processor_name_len);
    //MPI_Barrier(MPI_COMM_WORLD);
    printf("Worker %d of %d has started on %s with %d cpu cores.\n", workerId, world_size, processor_name, coresPerTask);
	fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    //Initial request to manager for assigning work item.
    MPI_Send(NULL, 0, MPI_INT, 0, WORK_REQ_TAG, MPI_COMM_WORLD);
    // Probe for an incoming message from process zero
    MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &wStatus);
    MPI_Get_count(&wStatus, MPI_CHAR, &nameLen);
    MPI_Recv(thisLigandName, nameLen, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &wStatus);
    thisLigandName[nameLen] = '\0';  // i don't know why, but on tianhe3, there is a 0xFF char at the end of thisLigandName. this line should fix that annoying problem
    //printf("INFO worker %d receive %d char: %s.\n", workerId, nameLen, thisLigandName);
    //printf("INFO worker %s receive char.\n", workerId);

    while (wStatus.MPI_TAG == COMPUTE_TAG)
    {
	// prepare ligand
        char fname_ligand[100]; //  = "/dev/shm/";
        sprintf(fname_ligand, "%s/worker_%d.in_ligand.sdf", TMPDIR, workerId);
        char sqlcmd[500];
        //strcpy(thisLigandName, "00003");
        sprintf(sqlcmd, "sqlite3 protn_lignd.db 'SELECT molstr FROM ligand WHERE molid = %s ;' > %s", thisLigandName, fname_ligand);
        //sprintf(sqlcmd, "sqlite3 protn_lignd.db 'SELECT molstr FROM ligand WHERE molid = 00003 ;' > %s", fname_ligand);
        //sprintf(sqlcmd, "time sqlite3 protn_lignd.db 'SELECT molstr FROM ligand WHERE molid = %s ;' > /dev/null", thisLigandName);
        printf("Worker = %d : sqlcmd = %s\n", workerId, sqlcmd);
        system(sqlcmd);
        sqlcmd[0] = '\0';

        char fname_protein[100]; // = "/dev/shm/";
        sprintf(fname_protein, "%s/worker_%d.in_protein.pdb", TMPDIR, workerId);
        sprintf(sqlcmd, "sqlite3 protn_lignd.db 'SELECT molstr FROM protein WHERE molid = %d ;' > %s", workerId % 99 + 1, fname_protein);
        //sprintf(sqlcmd, "time sqlite3 protn_lignd.db 'SELECT molstr FROM protein WHERE molid = %d ;' > /dev/null", workerId % 99 + 1);
        printf("Worker = %d : sqlcmd = %s\n", workerId, sqlcmd);
        system(sqlcmd);
        sqlcmd[0] = '\0';

        printf("Worker = %d : processing ligand '%s' ...\n", workerId, thisLigandName);
		fflush(stdout);

        char sminaCmd[500]; // TODO: how to detect number of cores assigned to each task? programmatically adjust --cpu
        sprintf(sminaCmd, "smina --seed 0 --cpu %d --num_modes 1 --receptor %s --ligand %s --autobox_ligand %s --out %s/worker_%d.out_ligand_%s.sdf --quiet || echo WARNING: smina not_ok", 
                          coresPerTask, fname_protein, fname_ligand, fname_ligand, TMPDIR, workerId, thisLigandName);
        // sprintf(sminaCmd, "time head -c 400M /dev/urandom | gzip > /dev/null");
        // sprintf(sminaCmd, "time gzip -c urandom1500M > /dev/null");
        //sprintf(sminaCmd, "time smina --seed 0 --cpu 1 --num_modes 1 --receptor /path/to/test_smina/5aba/5aba_protein.pdb --ligand /path/to/test_smina/5aba/5aba_ligand.sdf --autobox_ligand /path/to/test_smina/5aba/5aba_ligand.sdf --out /dev/shm/%d.sdf --quiet", workerId);
        //Ask smina to perform molecular docking.
        printf("Worker = %d : ligand '%s' sminaCmd=%s\n", workerId, thisLigandName, sminaCmd);
		//fflush(stdout);
        system(sminaCmd);

        sminaCmd[0] = '\0';

	// int a = 2;
	// float x = (float)rand()/(float)(RAND_MAX/a);
        // sleep(x);

        numJobs++;

        //Request for another work item.
        MPI_Send(NULL, 0, MPI_INT, 0, WORK_REQ_TAG, MPI_COMM_WORLD);
        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &wStatus);
        MPI_Get_count(&wStatus, MPI_CHAR, &nameLen);
        MPI_Recv(thisLigandName, nameLen, MPI_CHAR, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &wStatus);
        //printf("INFO worker %d receive %d char: %s.\n", workerId, nameLen, thisLigandName);
    }

    char tarcmd[500];
    sprintf(tarcmd, "tar cvzf ./output/worker_%d.tar.gz %s/worker_%d.out*.sdf && rm -rf %s/worker_%d.* ", workerId, TMPDIR, workerId, TMPDIR, workerId);
    printf("Worker = %d : archiving %s\n", workerId, tarcmd);
    //fflush(stdout);
    system(tarcmd);

    endTime = MPI_Wtime(); // end timer.
    if (wStatus.MPI_TAG == TERMINATE_TAG)
    {
        printf("Worker %d sminaMpiWorker exiting after finishing %d jobs in %.2lf seconds.\n", workerId, numJobs, endTime - startTime);
		fflush(stdout);
    }
    else
    {
        printf("Worker %d has received invalid Tag\n", workerId);
		fflush(stdout);
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    return;
}
