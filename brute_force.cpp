#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <eigen3/Eigen/Geometry>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <omp.h>
#include <unistd.h>
#include <vector>
using namespace std;

// additional parameters
double sphereRadius = 8;
int SPHERES;
int ATOMS;
int FRAMES;
int FRAMEONE;
int FRAMETWO;
// Atoms[<frame>][<atom>][<coordinate>]
vector<vector<vector<double>>> A;
// Maps sphere to CA; CAAtomNumber[<sphere>]
vector<int> sphereCA;
// List of atoms in [<sphere>]
vector<vector<int>>* sphereAtoms;

struct {
    int i = -1;
    int j = -1;
    double rmsd = -1;
} currentBest;

int AllocationsCountGlobal = 0;
int RMSDCalculationCountGlobal = 0;
int omp_thread_id;
omp_lock_t matrixMutex;

struct {
    string trajectory = "";
    int matrix_size = -1;
    int omp_threads_number = 1;
    bool generate_image = false;
} config;
//--------------------------------------------------

#pragma omp threadprivate(omp_thread_id, FRAMEONE, FRAMETWO)

typedef std::chrono::time_point<std::chrono::steady_clock, std::chrono::nanoseconds> timestamp;

class ProgressBar {
  public:
    ProgressBar(int totalSteps) : totalSteps(totalSteps), startTime(std::chrono::steady_clock::now()) {}

    void update(int currentStep) {
        int progress = (currentStep * 100) / totalSteps;
        if (progress == 0)
            return;
        auto currentTime = std::chrono::steady_clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count();
        int estimatedTotalTime = (elapsedTime * 100) / progress;
        int estimatedRemainingTime = estimatedTotalTime - elapsedTime;

        std::cout << "[" << std::string(progress / 5,  '#') << std::string((100 - progress) / 5, ' ') << "] " << std::fixed << std::setprecision(2)
                  << static_cast<double>(progress) << "%"
                  << " Time left: " << estimatedRemainingTime << "s     \r";
        std::cout.flush();
    }

  private:
    int totalSteps;
    std::chrono::steady_clock::time_point startTime;
};

class Progress {
  private:
    int barWidth;
    int currentSteps;
    int allStepsCount;

    // determine time left helper variables
    timestamp lastTimestamp;
    int timeLeftH;
    int timeLeftMin;
    int timeLeftSec;

    int modulo;

  public:
    Progress(int stepsCount, int mod = 100000) {
        allStepsCount = stepsCount;
        barWidth = 70;
        currentSteps = 0;

        lastTimestamp = std::chrono::steady_clock::now();

        modulo = mod;

        std::cout << "[>";
        for (int i = 1; i < barWidth; ++i) {
            std::cout << " ";
        }
        printf("] 0 %%\r");
        std::cout.flush();
    }

    void improve() {

        if (currentSteps % modulo) {
            currentSteps++;
            return;
        }
        std::cout << "[";
        int pos = barWidth * float(currentSteps) / float(allStepsCount);
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos)
                std::cout << "=";
            else if (i == pos)
                std::cout << ">";
            else
                std::cout << " ";
        }

        timestamp currentTimestamp = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(currentTimestamp - lastTimestamp);
        if (elapsed.count() > 5) {
            lastTimestamp = currentTimestamp;
            int stepsLeft = allStepsCount - currentSteps;

            timeLeftH = stepsLeft / modulo * elapsed / 1h;
            timeLeftMin = stepsLeft / modulo * elapsed / 1min;
            timeLeftSec = stepsLeft / modulo * elapsed / 1s;
        }
        printf("] %.2f%%  %dh:%dm:%ds left   \r", float(currentSteps) / float(allStepsCount) * 100.0, timeLeftH, timeLeftMin % 60, timeLeftSec % 60);
        std::cout.flush();
        currentSteps++;
    }

    void end() {
        std::cout << "[";
        int pos = barWidth * 1;
        for (int i = 0; i < barWidth; ++i) {
            if (i < pos)
                std::cout << "=";
            else if (i == pos)
                std::cout << ">";
            else
                std::cout << " ";
        }
        std::cout << "] 100%                       ";
        std::cout.flush();
        std::cout << std::endl;
        std::cout << "Done!" << std::endl;
    }
};

// reading data from input pdb file.
void readFile(string filename) {
    cout << "Reading file: " << filename << endl;
    string line;
    ifstream myfile1(filename);
    int lines_count = 0;
    if (myfile1.is_open()) {
        while (getline(myfile1, line)) {
            lines_count++;
        }
        myfile1.close();
    } else {
        cout << "Cannot find file: " << filename << endl;
        return;
    }
    Progress p1(lines_count);
    ifstream myfile(filename);
    if (myfile.is_open()) {
        int frame = 0;
        int atom;
        SPHERES = 0;
        FRAMES = 0;
        ATOMS = 0;
        A = {};
        sphereCA = {};
        while (getline(myfile, line)) {
            p1.improve();
            if (line[0] == 'M') {
                frame = stoi(line.substr(9, 5));
                frame--;
                A.push_back({});
                FRAMES++;
            } else if (line[0] == 'A') {
                atom = stoi(line.substr(6, 5));
                atom--;
                A[frame].push_back({});
                A[frame][atom].push_back(stod(line.substr(30, 8)));
                A[frame][atom].push_back(stod(line.substr(38, 8)));
                A[frame][atom].push_back(stod(line.substr(46, 8)));
                if (frame == 0) {
                    ATOMS++;
                    if (line[14] == 'A' and line[13] == 'C') {
                        sphereCA.push_back(atom);
                        SPHERES++;
                    }
                }
            }
        }
        p1.end();
        myfile.close();
        cout << "File parsed" << endl;
    } else {
        cout << "Cannot find file: " << filename << endl;
    }
}

// Find3DAffineTransform is from oleg-alexandrov repository on github, available here https://github.com/oleg-alexandrov/projects/blob/master/eigen/Kabsch.cpp
// [as of 27.01.2022] Given two sets of 3D points, find the rotation + translation + scale which best maps the first set to the second. Source:
// http://en.wikipedia.org/wiki/Kabsch_algorithm

// The input 3D points are stored as columns.
Eigen::Affine3d Find3DAffineTransform(Eigen::Matrix3Xd in, Eigen::Matrix3Xd out) {

    // Default output
    Eigen::Affine3d A;
    A.linear() = Eigen::Matrix3d::Identity(3, 3);
    A.translation() = Eigen::Vector3d::Zero();

    if (in.cols() != out.cols())
        throw "Find3DAffineTransform(): input data mis-match";

    // First find the scale, by finding the ratio of sums of some distances,
    // then bring the datasets to the same scale.
    double dist_in = 0, dist_out = 0;
    for (int col = 0; col < in.cols() - 1; col++) {
        dist_in += (in.col(col + 1) - in.col(col)).norm();
        dist_out += (out.col(col + 1) - out.col(col)).norm();
    }
    if (dist_in <= 0 || dist_out <= 0)
        return A;
    double scale = dist_out / dist_in;
    out /= scale;

    // Find the centroids then shift to the origin
    Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
    Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
    for (int col = 0; col < in.cols(); col++) {
        in_ctr += in.col(col);
        out_ctr += out.col(col);
    }
    in_ctr /= in.cols();
    out_ctr /= out.cols();
    for (int col = 0; col < in.cols(); col++) {
        in.col(col) -= in_ctr;
        out.col(col) -= out_ctr;
    }

    // SVD
    Eigen::MatrixXd Cov = in * out.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Find the rotation
    double d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    if (d > 0)
        d = 1.0;
    else
        d = -1.0;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
    I(2, 2) = d;
    Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();

    // The final transform
    A.linear() = scale * R;
    A.translation() = scale * (out_ctr - R * in_ctr);

    return A;
}

// superpose changes vector of atoms frame 2 to map atoms from frame 1 in the way to minimise RMSD between both frames
void superpose(const vector<vector<double>>& frame1, vector<vector<double>>& frame2) {
    int atomsInSphere = frame1.size();
    Eigen::Matrix3Xd S1(3, atomsInSphere);
    Eigen::Matrix3Xd S2(3, atomsInSphere);
    for (int j = 0; j < atomsInSphere; j++) {
        for (int k = 0; k < 3; k++) {
            S1(k, j) = frame1[j][k];
            S2(k, j) = frame2[j][k];
        }
    }
    Eigen::Affine3d RT = Find3DAffineTransform(S2, S1);
    S2 = RT.linear() * S2;
    for (int j = 0; j < atomsInSphere; j++) {
        S2.block<3, 1>(0, j) += RT.translation();
    }
    for (int j = 0; j < atomsInSphere; j++) {
        for (int k = 0; k < 3; k++) {
            frame2[j][k] = S2(k, j);
        }
    }
}

// calculating distance between 2 atoms, used when allocating atoms to spheres.
double atomsDistanceCalc(int atom1, int atom2) {
    double result = (A[FRAMEONE][atom1][0] - A[FRAMEONE][atom2][0]) * (A[FRAMEONE][atom1][0] - A[FRAMEONE][atom2][0]) +
                    (A[FRAMEONE][atom1][1] - A[FRAMEONE][atom2][1]) * (A[FRAMEONE][atom1][1] - A[FRAMEONE][atom2][1]) +
                    (A[FRAMEONE][atom1][2] - A[FRAMEONE][atom2][2]) * (A[FRAMEONE][atom1][2] - A[FRAMEONE][atom2][2]);
    return sqrt(result);
}

// allocating atoms into spheres, based on sphereRadius
void atomsAllocation(int firstFrame) {
    FRAMEONE = firstFrame;
    std::vector<int> temp;
    temp.reserve(ATOMS);
    sphereAtoms[omp_thread_id].assign(SPHERES, temp);
    for (int i = 0; i < ATOMS; i++) {
        for (int j = 0; j < SPHERES; j++) {
            if (atomsDistanceCalc(i, sphereCA[j]) <= sphereRadius) {
                sphereAtoms[omp_thread_id][j].push_back(i);
            }
        }
    }
}

// calculating RMSD on spheres, on choosen frames
double calculateRMSDSuperpose(int secondFrame) {
    FRAMETWO = secondFrame;
    double result = 0;
    std::vector<std::vector<std::vector<double>>> sphereMatrix;
    for (int s = 0; s < SPHERES; s++) {
        int atomsInSphere = sphereAtoms[omp_thread_id][s].size();
        sphereMatrix.assign(2, {});
        sphereMatrix[0].assign(atomsInSphere, {});
        sphereMatrix[1].assign(atomsInSphere, {});

        for (int j = 0; j < atomsInSphere; j++) {
            sphereMatrix[0][j] = A[FRAMEONE][sphereAtoms[omp_thread_id][s][j]];
            sphereMatrix[1][j] = A[FRAMETWO][sphereAtoms[omp_thread_id][s][j]];
        }
        double tempResult = 0;
        superpose(sphereMatrix[0], sphereMatrix[1]);
        for (int j = 0; j < atomsInSphere; j++) {
            for (int k = 0; k < 3; k++) {
                double tempRMSD = sphereMatrix[1][j][k] - sphereMatrix[0][j][k];
                tempResult += tempRMSD * tempRMSD;
            }
        }
        tempResult /= atomsInSphere * 3.0;
        tempResult = sqrt(tempResult);
        result += tempResult;
    }
    return result;
}

string getCurrentDateTime() {
    auto now = chrono::system_clock::now();
    time_t currentTime = chrono::system_clock::to_time_t(now);
    stringstream ss;
    ss << put_time(localtime(&currentTime), "%d-%m-%Y_%H-%M-%S");
    return ss.str();
}

string transformFilename(const string& filePath) {

    filesystem::path file_path(filePath);

    // Get the filename without the extension
    string filename = file_path.stem().string();
    string ext = file_path.extension().string();

    // Find the last dot in the filename
    size_t lastDotPos = ext.find_last_of('.');
    if (lastDotPos != string::npos) {
        // Replace the last dot with an underscore
        ext.replace(lastDotPos, 1, "_");
    }

    return filename + ext;
}

int main(int argc, char* argv[]) {
    omp_init_lock(&matrixMutex);
    if (argc == 1) {
        cout << "Too few arguments" << endl;
        cout << "Usage:" << endl;
        cout << "  brute_force <TRAJECTORY> [OPTIONS]..." << endl;
        cout << endl;
        cout << "Options:" << endl;
        cout << "  -s SIZE              setting matrix size to SIZE, default is max size for current trajectory" << endl;
        cout << "  -t THREADS           setting number of omp threads to THREADS, default is 1" << endl;
        cout << "  -i                   generating matrix image in gray scale:" << endl;
        cout << "                        - max_rmsd: white pixel," << endl;
        cout << "                        - min_rmsd: black pixel," << endl;
        cout << "                        - values are normalized." << endl;
        cout << "                        By default no image generated" << endl;
        cout << endl;
        cout << "Examples:" << endl;
        cout << "  brute_force ./trajectory.pdb               max matrix size and 1 omp thread" << endl;
        cout << "  brute_force ./trajectory.pdb -t 10 -i      using 10 omp threads and generating an image" << endl;
        return 0;
    } else {
        int option;

        while ((option = getopt(argc, argv, "s:t:i")) != -1) {
            switch (option) {
            case 's':
                config.matrix_size = std::stoi(optarg);
                break;
            case 't':
                config.omp_threads_number = std::stoi(optarg);
                break;
            case 'i':
                config.generate_image = true;
                break;
            case '?':
                // Handle invalid options or missing argument error
                std::cerr << "Invalid option: -" << optopt << std::endl;
                return 1;
            default:
                break;
            }
        }
        if (optind < argc) {
            config.trajectory = argv[optind];
        } else {
            // If the filename is missing, print an error and exit
            std::cerr << "Error: Missing filename!" << std::endl;
            return 1;
        }
    }

    cout << "Config:" << endl;
    cout << "config.trajectory: " << config.trajectory << endl;
    cout << "config.matrix_size: " << config.matrix_size << endl;
    cout << "config.omp_threads_number: " << config.omp_threads_number << endl;
    cout << "config.generate_image: " << boolalpha << config.generate_image << endl;

    readFile(config.trajectory);

    if (config.matrix_size == -1) {
        config.matrix_size = FRAMES;
    }
    cout << "updated: config.matrix_size: " << config.matrix_size << endl;
    int size = config.matrix_size;

    double** matrix = new double*[size];
    for (int i = 0; i < size; i++) {
        matrix[i] = new double[size];
    }

    omp_set_num_threads(config.omp_threads_number);

    cout << "[OMP] [Number of threads]: " << config.omp_threads_number << endl;
    sphereAtoms = new std::vector<std::vector<int>>[config.omp_threads_number];
    cout << "Calculating.." << endl;

    Progress p2(size * size / config.omp_threads_number, 100);
    auto start = std::chrono::steady_clock::now();

#pragma omp parallel for reduction(+ : AllocationsCountGlobal, RMSDCalculationCountGlobal)
    for (int i = 0; i < size; i++) {
        omp_thread_id = omp_get_thread_num();
        atomsAllocation(i);
        AllocationsCountGlobal++;
        int RMSDCalculationCount = 0;
        for (int j = 0; j < size; j++) {
            double result = calculateRMSDSuperpose(j);
            matrix[i][j] = result;
#pragma omp critical
            {
                if (i != j && result > currentBest.rmsd) {
                    cout << "[New best] [" << i << ", " << j << "] = " << result << endl;
                    currentBest = { i, j, result };
                }
            }
            if (omp_get_thread_num() == 0) {
                p2.improve();
            }
            RMSDCalculationCount++;
        }
        RMSDCalculationCountGlobal += RMSDCalculationCount;
    }

    p2.end();

    double max_value = matrix[0][0];
    int max_i = 0;
    int max_j = 0;

    // Find the maximum value in the matrix
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (matrix[i][j] > max_value) {
                max_value = matrix[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }
    auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    int hours = std::chrono::duration_cast<std::chrono::hours>(elapsed).count();
    int minutes = std::chrono::duration_cast<std::chrono::minutes>(elapsed).count() % 60;
    int seconds = std::chrono::duration_cast<std::chrono::seconds>(elapsed).count() % 60;

    cout << "Max value: " << max_value << endl;
    cout << "(i, j): (" << max_i << ", " << max_j << ")" << endl;
    cout << "Total rmsd calculated: " << RMSDCalculationCountGlobal << endl;
    cout << "Total allocations: " << AllocationsCountGlobal << endl;
    cout << "The calculation took: ";
    if (hours) {
        cout << hours << "h ";
    }
    if (minutes) {
        cout << minutes << "m ";
    }
    cout << seconds << "s" << endl;

    // creating an image

    if (config.generate_image) {

        // normalize the matrix
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                matrix[i][j] /= max_value;
            }
        }

        const int max_intensity = 255;

        string map_filename = "./" + transformFilename(config.trajectory) + "_" + getCurrentDateTime() + "_map.pgm";

        std::ofstream image_file(map_filename, std::ios::binary);

        // Write the PGM header
        image_file << "P5\n" << size << " " << size << "\n" << max_intensity << "\n";

        // Write the pixel values
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                unsigned char intensity = static_cast<unsigned char>((1 - matrix[i][j]) * max_intensity);
                image_file.write(reinterpret_cast<char*>(&intensity), sizeof(unsigned char));
            }
        }

        // Close the image file
        image_file.close();

        cout << "Map generated under: " << map_filename << endl;
    }
    for (int i = 0; i < size; i++) {
        delete[] matrix[i];
    }
    delete[] matrix;
    delete[] sphereAtoms;
    omp_destroy_lock(&matrixMutex);
    return 0;
}
