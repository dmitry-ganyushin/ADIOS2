//
// Created by ganyush on 9/29/20.
//

#include <adios2.h>
#include <chrono>
#include <ios>      //std::ios_base::failure
#include <iostream> //std::cout
#include <queue>
#include <set>
#include <stack>
#include <stdexcept> //std::invalid_argument std::exception
#include <vector>

#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

bool endsWith(const std::string &s, const std::string &suffix)
{
    return s.find(suffix, s.size() - suffix.size()) != std::string::npos;
}

bool startsWith(const std::string &s, const std::string &prefix)
{
    return s.rfind(prefix, 0) == 0;
}

void TraverseBreadthSearchFirst(adios2::Group &g)
{
    int numberOfGroups = 0;
    int numberOfVariables = 0;
    std::queue<std::string> queue;
    std::string root = g.InquirePath();
    queue.push(root);

    while (!queue.empty())
    {
        int levelSize = queue.size();
        for (int i = 0; i < levelSize; i++)
        {
            std::string currentNode = queue.front();
            queue.pop();
            // add the node to the current level
            g.setPath(currentNode);
            std::vector<std::string> groups = g.AvailableGroups();

            numberOfGroups += groups.size();

            std::vector<std::string> variables = g.AvailableVariables();
            numberOfVariables += variables.size();

            // insert the children of current node in the queue
            for (auto v : groups)
            {
                queue.push(currentNode + "/" + v);
            }
        }
    }
    std::cout << "Number of groups " << numberOfGroups << std::endl;
    std::cout << "Number of variables " << numberOfVariables << std::endl;
    return;
}

void TraverseDepthSearchFirst(adios2::Group &g)
{
    int numberOfGroups = 0;
    int numberOfVariables = 0;
    std::stack<std::string> stack;
    std::string root = g.InquirePath();
    stack.push(root);

    while (!stack.empty())
    {

        std::set<std::string> visited;

        // Pop a vertex from stack and print it
        std::string currentNode = stack.top();
        stack.pop();

        if (visited.find(currentNode) == visited.end())
        {
            g.setPath(currentNode);
            std::vector<std::string> groups = g.AvailableGroups();
            numberOfGroups += groups.size();

            std::vector<std::string> variables = g.AvailableVariables();
            numberOfVariables += variables.size();
            visited.insert(currentNode);

            // insert the children of current node in the queue
            for (auto v : groups)
            {
                stack.push(currentNode + "/" + v);
            }
        }
    }
    std::cout << "Number of groups " << numberOfGroups << std::endl;
    std::cout << "Number of variables " << numberOfVariables << std::endl;
    return;
}

int main(int argc, char *argv[])
{
    int rank, size;
#if ADIOS2_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#else
    rank = 0;
    size = 1;
#endif
    const int number_of_steps = 5;
    /** Application variable */
    std::string filename = "myVector_cppz4.bp";
    std::vector<float> myFloats = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<int> myInts = {0, -1, -2, -3, -4, -5, -6, -7, -8, -9};
    const std::size_t Nx = myFloats.size();

    const std::string myString("Hello Variable String from rank " +
                               std::to_string(rank));

    try
    {
        /** ADIOS class factory of IO class objects */
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif

        /*** IO class object: settings and factory of Settings: Variables,
         * Parameters, Transports, and Execution: Engines */
        adios2::IO bpIO = adios.DeclareIO("BPFile_N2N");
        /** Engine derived class, spawned to start IO operations */
        adios2::Engine bpFileWriter = bpIO.Open(filename, adios2::Mode::Append);

        /** global array : name, { shape (total) }, { start (local) }, {
         * count
         * (local) }, all are constant dimensions */

        adios2::Variable<float> bpFloats =
            bpIO.DefineVariable<float>("group1/subgroup1/bpFloats/__data__", {size * Nx},
                                       {rank * Nx}, {Nx}, adios2::ConstantDims);
        adios2::Variable<int> bpInts =
            bpIO.DefineVariable<int>("group1/subgroup1/bpInts/__data__", {size * Nx},
                                     {rank * Nx}, {Nx}, adios2::ConstantDims);

        for (int i = 0; i < number_of_steps; i++)
        {
            bpFileWriter.BeginStep();

            /** Put variables for buffering, template type is optional */
            if (i%2 == 0){
                bpFileWriter.Put<float>(bpFloats, myFloats.data());
            }


            bpFileWriter.Put(bpInts, myInts.data());

            bpFileWriter.EndStep();
        }
        /** Create bp file, engine becomes unreachable after this*/
        bpFileWriter.Close();
        if (rank == 0)
        {
            std::cout << "Wrote file " << filename << " to disk.  "
                      << std::endl;
        }
    }
    catch (std::invalid_argument &e)
    {
        std::cerr << "Invalid argument exception: " << e.what() << "\n";
#if ADIOS2_USE_MPI
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    catch (std::ios_base::failure &e)
    {
        std::cerr << "IO System base failure exception: " << e.what() << "\n";
#if ADIOS2_USE_MPI
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    catch (std::exception &e)
    {
        std::cerr << "Exception: " << e.what() << "\n";
#if ADIOS2_USE_MPI
        std::cerr << "STOPPING PROGRAM from rank " << rank << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }

    try
    {
        /** ADIOS class factory of IO class objects */
        // adios2::ADIOS adios(MPI_COMM_WORLD);
#if ADIOS2_USE_MPI
        adios2::ADIOS adios(MPI_COMM_WORLD);
#else
        adios2::ADIOS adios;
#endif

        /*** IO class object: settings and factory of Settings: Variables,
         * Parameters, Transports, and Execution: Engines */
        adios2::IO bpIO = adios.DeclareIO("ReadBP");

        const auto attributesInfo = bpIO.AvailableAttributes();

        for (const auto &attributeInfoPair : attributesInfo)
        {
            std::cout << "Attribute: " << attributeInfoPair.first;
            for (const auto &attributePair : attributeInfoPair.second)
            {
                std::cout << "\tKey: " << attributePair.first
                          << "\tValue: " << attributePair.second << "\n";
            }
            std::cout << "\n";
        }

        adios2::Engine bpReader = bpIO.Open(filename, adios2::Mode::Read);

        for (int step = 0; step < number_of_steps; step++)
        {
            bpReader.BeginStep();

            std::cout << "****Step " << step << "*** " << std::endl;
            adios2::Group g = bpIO.GetGroup("group1");
            std::cout << "****Path of g***" << g.InquirePath() << std::endl;
            std::cout << "****Available variables*********" << std::endl;
            auto result = g.AvailableVariables();
            for (auto v : result)
                std::cout << v << " ";
            std::cout << std::endl
                      << "********************************" << std::endl;

            std::cout << "****Available groups************" << std::endl;
            result = g.AvailableGroups();
            for (auto v : result)
                std::cout << v << " ";
            std::cout << std::endl
                      << "********************************" << std::endl;

            adios2::Group g1 = g.OpenGroup("subgroup1");
            std::cout << "****Path of g1***" << g1.InquirePath() << std::endl;
            std::cout << "****Available variables*********" << std::endl;
            result = g1.AvailableVariables();
            for (auto v : result)
                std::cout << v << " ";
            std::cout << std::endl
                      << "********************************" << std::endl;

            std::cout << "****Available groups************" << std::endl;
            result = g1.AvailableGroups();
            for (auto v : result)
                std::cout << v << " ";
            std::cout << std::endl
                      << "********************************" << std::endl;


            adios2::Group g2 = g1.OpenGroup("bpFloats");
            std::cout << "****Path of g2 ***" << g2.InquirePath() << std::endl;
            std::cout << "****Available variables*********" << std::endl;
            result = g2.AvailableVariables();
            for (auto v : result)
                std::cout << v << " ";
            std::cout << std::endl
                      << "********************************" << std::endl;

            std::cout << "****Available groups************" << std::endl;
            result = g2.AvailableGroups();
            for (auto v : result)
                std::cout << v << " ";
            std::cout << std::endl
                      << "********************************" << std::endl;


            adios2::Group g3 = g1.OpenGroup("bpInts");
            std::cout << "****Path of g3 ***" << g3.InquirePath() << std::endl;
            std::cout << "****Available variables*********" << std::endl;
            result = g3.AvailableVariables();
            for (auto v : result)
                std::cout << v << " ";
            std::cout << std::endl
                      << "********************************" << std::endl;

            std::cout << "****Available groups************" << std::endl;
            result = g3.AvailableGroups();
            for (auto v : result)
                std::cout << v << " ";
            std::cout << std::endl
                      << "********************************" << std::endl;

            /** Write variable for buffering */
            adios2::Variable<float> bpFloats =
                g2.InquireVariable<float>("__data__");
            adios2::Variable<int> bpInts = g3.InquireVariable<int>("__data__");

            const std::size_t Nx = 10;
            if (bpFloats) // means found
            {
                std::vector<float> myFloats;

                // read only the chunk corresponding to our rank
                bpFloats.SetSelection({{Nx * rank}, {Nx}});
                // myFloats.data is pre-allocated
                bpReader.Get<float>(bpFloats, myFloats, adios2::Mode::Sync);

                if (rank == 0)
                {
                    std::cout << "MyFloats: \n";
                    for (const auto number : myFloats)
                    {
                        std::cout << number << " ";
                    }
                    std::cout << "\n";
                }
            }

            if (bpInts) // means not found
            {
                std::vector<int> myInts;
                // read only the chunk corresponding to our rank
                bpInts.SetSelection({{Nx * rank}, {Nx}});

                bpReader.Get<int>(bpInts, myInts, adios2::Mode::Sync);

                if (rank == 0)
                {
                    std::cout << "myInts: \n";
                    for (const auto number : myInts)
                    {
                        std::cout << number << " ";
                    }
                    std::cout << "\n";
                }
            }

            bpReader.EndStep();
        }
        /** Close bp file, engine becomes unreachable after this*/
        bpReader.Close();
    }
    catch (std::invalid_argument &e)
    {
        if (rank == 0)
        {
            std::cerr
                << "Invalid argument exception, STOPPING PROGRAM from rank "
                << rank << "\n";
            std::cerr << e.what() << "\n";
        }
#if ADIOS2_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    catch (std::ios_base::failure &e)
    {
        if (rank == 0)
        {
            std::cerr << "IO System base failure exception, STOPPING PROGRAM "
                         "from rank "
                      << rank << "\n";
            std::cerr << e.what() << "\n";
            std::cerr << "The file " << filename << " does not exist."
                      << " Presumably this is because hello_bpWriter has not "
                         "been run."
                      << " Run ./hello_bpWriter before running this program.\n";
        }
#if ADIOS2_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
    catch (std::exception &e)
    {
        if (rank == 0)
        {
            std::cerr << "Exception, STOPPING PROGRAM from rank " << rank
                      << "\n";
            std::cerr << e.what() << "\n";
        }
#if ADIOS2_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }

#if ADIOS2_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
