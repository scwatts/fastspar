#include "bootstrap_opts.h"


void print_help() {
    fprintf(stderr, "Program: FastSpar bootstraps (c++ implementation of SparCC)\n");
    fprintf(stderr, "Version 0.0.4\n");
    fprintf(stderr, "Contact: Stephen Watts (s.watts2@student.unimelb.edu.au)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  bootstraps --otu_table <of> --number <n> --prefix <p>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -c/--otu_table <of>   OTU input table\n");
    fprintf(stderr, "  -n/--number <n>       Number of bootstrap samples to generate\n");
    fprintf(stderr, "  -p/--prefix <p>       Prefix out bootstrap output files\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -t <th>, --threads <th>\n");
    fprintf(stderr, "                        Number of threads (1 default)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Other:\n");
    fprintf(stderr, "  -h        --help\n");
    fprintf(stderr, "                Display this help and exit\n");
    fprintf(stderr, "  -v        --version\n");
    fprintf(stderr, "                Display version information and exit\n");
}


void print_version() {
    fprintf(stderr, "Program: FastSpar bootstraps (c++ implementation of SparCC)\n");
    fprintf(stderr, "Version: 0.0.4\n");
    fprintf(stderr, "Contact: Stephen Watts (s.watts2@student.unimelb.edu.au)\n");
}


BootstrapOptions get_commandline_arguments(int argc, char **argv) {
    // Get instance of BootstrapOptions
    BootstrapOptions bootstrap_options;

    // Commandline arguments (for getlongtops)
    struct option long_options[] =
        {
            {"otu_table", required_argument, NULL, 'c'},
            {"number", required_argument, NULL, 'n'},
            {"threads", required_argument, NULL, 't'},
            {"prefix", required_argument, NULL, 'p'},
            {"version", no_argument, NULL, 'v'},
            {"help", no_argument, NULL, 'h'},
            {NULL, 0, 0, 0}
        };

    // Parse commandline arguments
    while (1) {
        // Parser variables
        int option_index = 0;
        int c;

        // Parser
        c = getopt_long(argc, argv, "hvc:p:n:t:", long_options, &option_index);

        // If no more arguments to parse, break
        if (c == -1) {
            break;
        }

        switch(c) {
            case 'c':
                bootstrap_options.otu_filename = optarg;
                break;
            case 'p':
                bootstrap_options.bootstrap_prefix = optarg;
                break;
            case 'n':
                bootstrap_options.bootstrap_number = int_from_optarg(optarg);
                break;
            case 't':
                bootstrap_options.threads = int_from_optarg(optarg);
                break;
            case 'v':
                print_version();
                exit(0);
            case 'h':
                print_help();
                exit(0);
            default:
                exit(1);
        }
    }


    // Check if have an attempt at arguments
    if (argc < 7) {
        print_help();
        fprintf(stderr,"\n%s: error: option -c/--otu_table, -n/--number, and -p/--prefix are required\n", argv[0]);
        exit(1);
    }

    // Abort execution if given unknown arguments
    if (optind < argc){
        print_help();
        fprintf(stderr, "\n%s: invalid argument: %s\n", argv[0], argv[optind++]);
    }


    // Make sure we have filenames and parameters
    if (bootstrap_options.otu_filename.empty()) {
        print_help();
        fprintf(stderr,"\n%s: error: argument -c/--otu_table is required\n", argv[0]);
        exit(1);
    }
    if (bootstrap_options.bootstrap_prefix.empty()) {
        print_help();
        fprintf(stderr,"\n%s: error: argument -p/--prefix is required\n", argv[0]);
        exit(1);
    }
    if (bootstrap_options.bootstrap_number == 0) {
        print_help();
        fprintf(stderr,"\n%s: error: argument -n/--number is required\n", argv[0]);
        exit(1);
    }


    // Check if have a reasonable number of threads
    if (bootstrap_options.threads < 1) {
        print_help();
        fprintf(stderr,"\n%s: error: must specify at least 1 thread\n", argv[0]);
        exit(1);
    }

    // Check we don't attempt to use more threads than we have
    unsigned int available_threads = std::thread::hardware_concurrency();
    if (available_threads > 1 && bootstrap_options.threads > available_threads) {
        print_help();
        fprintf(stderr, "\n%s: error: only %d threads are available\n", argv[0], available_threads);
        exit(1);
    } else if (bootstrap_options.threads > 64) {
        print_help();
        fprintf(stderr, "\n%s: error: current hardcode limit of 64 threads\n", argv[0]);
        exit(1);
    }


    // Check that the OTU file exists
    std::ifstream otu_file;
    otu_file.open(bootstrap_options.otu_filename);
    if (!otu_file.good()) {
        print_help();
        fprintf(stderr, "\n%s: error: OTU table %s does not exist\n", argv[0], bootstrap_options.otu_filename.c_str());
        exit(1);
    }

    return bootstrap_options;
}
