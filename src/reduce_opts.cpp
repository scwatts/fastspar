#include "reduce_opts.h"


void print_help() {
    fprintf(stderr, "Program: Filter correlation and p-value table. Output as sparse matrices\n");
    fprintf(stderr, "Version 0.0.3\n");
    fprintf(stderr, "Contact: Stephen Watts (s.watts2@student.unimelb.edu.au)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:\n");
    fprintf(stderr, "  reduce --correlation_table <cf> --pvalue_table <pf> --correlation <rt> --pvalue <pt> --output_prefix <op>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  -r/--correlation_table <cf>       Correlation input table\n");
    fprintf(stderr, "  -p/--pvalue_table <pf>            Number of bootstrap samples to generate\n");
    fprintf(stderr, "  -o/--output_prefix <op>           Output prefix\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -t/--correlation <rt>             Absolute (sign is ignored) correlation threshold (default: 0.1)\n");
    fprintf(stderr, "  -l/--pvalue <pt>                  P-value threshold (default: 0.05)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Other:\n");
    fprintf(stderr, "  -h        --help\n");
    fprintf(stderr, "                Display this help and exit\n");
    fprintf(stderr, "  -v        --version\n");
    fprintf(stderr, "                Display version information and exit\n");
}


void print_version() {
    fprintf(stderr, "Program: Filter correlation and p-value table. Output as sparse matrices\n");
    fprintf(stderr, "Version: 0.0.3\n");
    fprintf(stderr, "Contact: Stephen Watts (s.watts2@student.unimelb.edu.au)\n");
}


ReduceOptions get_commandline_arguments(int argc, char **argv) {
    // Get instance of ReduceOptions
    ReduceOptions reduce_options;

    // Commandline arguments (for getlongtops)
    struct option long_options[] =
        {
            {"correlation_table", required_argument, NULL, 'r'},
            {"pvalue_table", required_argument, NULL, 'p'},
            {"correlation", required_argument, NULL, 'r'},
            {"pvalue", required_argument, NULL, 'l'},
            {"output_prefix", required_argument, NULL, 'o'},
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
        c = getopt_long(argc, argv, "hvc:p:r:l:o:", long_options, &option_index);

        // If no more arguments to parse, break
        if (c == -1) {
            break;
        }

        // Process current argument
        switch(c) {
            case 'r':
                reduce_options.correlation_filename = optarg;
                break;
            case 'p':
                reduce_options.pvalue_filename = optarg;
                break;
            case 't':
                reduce_options.correlation_threshold = float_from_optarg(optarg);
                break;
            case 'l':
                reduce_options.pvalue_threshold = float_from_optarg(optarg);
                break;
            case 'o':
                reduce_options.output_prefix = optarg;
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
        fprintf(stderr,"\n%s: error: option -r/--correlation_table, -p/--pvalue_table, and -o/--output_prefix are required\n", argv[0]);
        exit(1);
    }

    // Abort execution if given unknown arguments
    if (optind < argc){
        print_help();
        fprintf(stderr, "\n%s: invalid argument: %s\n", argv[0], argv[optind++]);
    }


    // Make sure we have filenames and parameters
    if (reduce_options.correlation_filename.empty()) {
        print_help();
        fprintf(stderr, "\n%s: error: argument -c/--correlation_table is required\n", argv[0]);
        exit(1);
    }
    if (reduce_options.pvalue_filename.empty()) {
        print_help();
        fprintf(stderr, "\n%s: error: argument -p/--pvalue_table is required\n", argv[0]);
        exit(1);
    }
    if (reduce_options.output_prefix.empty()) {
        print_help();
        fprintf(stderr, "\n%s: error: argument -o/--output_prefixis required\n", argv[0]);
        exit(1);
    }


    // Check that the correlation file exists
    std::ifstream checkfile;
    checkfile.open(reduce_options.correlation_filename);
    if (!checkfile.good()) {
        fprintf(stderr, "\n%s: error: correlation table file %s does not exist\n", argv[0], reduce_options.correlation_filename.c_str());
        exit(1);
    }
    checkfile.close();

    // Check that the pvalue file exists
    checkfile.open(reduce_options.pvalue_filename);
    if (!checkfile.good()) {
        fprintf(stderr, "\n%s: error: pvalue table file %s does not exist\n", argv[0], reduce_options.pvalue_filename.c_str());
        exit(1);
    }


    return reduce_options;
}
