#include <vector>
#include <string>

#define infinite 9999999
struct Parameters
{
    static Parameters& get(){
        static Parameters param;
        return param;
    }
    void set(size_t param, std::string param_read)
    {
        if (param_read == "KmerSize")
            Parameters::get().kmerSize = param;
        if (param_read == "Accumulative")
            Parameters::get().accumulative_h = param;
    }
    size_t accumulative_h;
    size_t kmerSize;
    size_t numThreads;
};

struct Progress
{
    static void update(size_t num_actual){
        if (Progress::get().show) {
            int val=(int) (((float) num_actual) / ((float) Progress::get().size_total) * 100);
            if (val % 10 == 0)
                if (!Progress::get().f[val/10]) {
                    Progress::get().f[val/10] = true;
                    show_progress(num_actual);
                }
        }
    }

    static Progress& get(){
        static Progress progress;
        return progress;
    }
    std::vector<bool> f = std::vector<bool>(10,false);
    size_t size_total;
    bool show = false;
private:
    static void show_progress(size_t num_actual){
        size_t val = (size_t)((float)num_actual/
                              (float)Progress::get().size_total * 100);
        std::cout << val << ((val==100)?"%\n":"% ") << std::flush;
        if (val == 100)
            _prepare_next();
    }
    static void _prepare_next(){
        Progress::get().f = std::vector<bool>(10,false);
    }
};
