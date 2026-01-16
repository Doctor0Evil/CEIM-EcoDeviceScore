#ifndef CEIM_ECODEVICESCORE_HPP
#define CEIM_ECODEVICESCORE_HPP

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <cctype>
#include <map>

namespace ceim
{
    // -----------------------------
    // Core CEIM node impact types
    // -----------------------------

    struct CeimNodeSample
    {
        // Concentrations (in consistent units, e.g., mg/L)
        double cin;   // Cin,x(t)
        double cout;  // Cout,x(t)

        // Flow Q(t) (e.g., m^3/s)
        double flow;

        // Time step length dt (seconds)
        double dt;

        // Regulatory / health reference concentration Cref,x
        double cref;

        // Hazard weight lambda_x (dimensionless)
        double lambda;
    };

    struct CeimNodeImpactResult
    {
        // Raw integral: sum_x lambda_x ∫ (Cin - Cout) * Q / Cref dt
        double kn;

        // Optional uncertainty; here kept simple as a placeholder
        double kn_std;
    };

    inline double compute_node_impact(
        const std::vector<CeimNodeSample>& samples
    )
    {
        double kn = 0.0;
        for (const auto& s : samples)
        {
            if (s.cref <= 0.0)
            {
                continue; // skip invalid reference
            }

            const double delta_c = s.cin - s.cout;  // Cin - Cout
            const double load_term = delta_c * s.flow * s.dt; // mass-like term

            const double normalized = (load_term / s.cref) * s.lambda;
            kn += normalized;
        }
        return kn;
    }

    inline CeimNodeImpactResult compute_node_impact_with_uncertainty(
        const std::vector<CeimNodeSample>& samples
    )
    {
        const double kn = compute_node_impact(samples);

        // Simple placeholder: no detailed propagation, but keep field present.
        CeimNodeImpactResult result;
        result.kn = kn;
        result.kn_std = 0.0;
        return result;
    }

    // ---------------------------------------------------
    // EcoDevice lifecycle score / CSV ingest
    // ---------------------------------------------------

    struct EcoDeviceRecord
    {
        std::string device_id;
        std::string category;
        bool        software_first;
        double      embodied_co2_kg;
        double      expected_lifetime_years;
        double      annual_energy_savings_kwh;
        double      annual_pollutant_reduction_tons;
        double      ecoimpact_score_01;
        std::string notes;
    };

    struct EcoDeviceCorridorAggregate
    {
        // Simple corridor/group ID key
        std::string corridor_id;

        // Aggregated metrics
        double total_embodied_co2_kg = 0.0;
        double total_annual_energy_savings_kwh = 0.0;
        double total_annual_pollutant_reduction_tons = 0.0;

        // Average ecoimpact score (0-1) across devices
        double average_ecoimpact_score_01 = 0.0;

        // Number of devices in corridor
        std::size_t device_count = 0;
    };

    namespace detail
    {
        inline std::string trim(const std::string& s)
        {
            std::size_t start = 0;
            while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start])))
            {
                ++start;
            }

            std::size_t end = s.size();
            while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1])))
            {
                --end;
            }

            return s.substr(start, end - start);
        }

        inline std::vector<std::string> split_csv_line(const std::string& line)
        {
            std::vector<std::string> result;
            std::string current;
            bool in_quotes = false;

            for (char ch : line)
            {
                if (ch == '"')
                {
                    in_quotes = !in_quotes;
                    continue;
                }

                if (ch == ',' && !in_quotes)
                {
                    result.push_back(trim(current));
                    current.clear();
                }
                else
                {
                    current.push_back(ch);
                }
            }

            result.push_back(trim(current));
            return result;
        }

        inline double parse_double(const std::string& s)
        {
            if (s.empty())
            {
                return 0.0;
            }
            char* end = nullptr;
            const double val = std::strtod(s.c_str(), &end);
            if (end == s.c_str())
            {
                return 0.0;
            }
            return val;
        }

        inline bool parse_bool(const std::string& s)
        {
            const std::string lower = [&]() {
                std::string tmp;
                tmp.reserve(s.size());
                for (char c : s)
                {
                    tmp.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
                }
                return tmp;
            }();

            return (lower == "1" || lower == "true" || lower == "yes");
        }
    } // namespace detail

    inline std::vector<EcoDeviceRecord> load_eco_device_csv(
        const std::string& filepath
    )
    {
        std::ifstream file(filepath);
        if (!file.is_open())
        {
            throw std::runtime_error("Failed to open CSV file: " + filepath);
        }

        std::string header_line;
        if (!std::getline(file, header_line))
        {
            throw std::runtime_error("CSV file is empty: " + filepath);
        }

        const auto header_cols = detail::split_csv_line(header_line);
        using IndexMap = std::map<std::string, std::size_t>;
        IndexMap idx;

        for (std::size_t i = 0; i < header_cols.size(); ++i)
        {
            idx[header_cols[i]] = i;
        }

        auto get_index = [&](const std::string& name) -> std::size_t {
            auto it = idx.find(name);
            if (it == idx.end())
            {
                throw std::runtime_error("Missing required column in CSV: " + name);
            }
            return it->second;
        };

        const std::size_t i_device_id                    = get_index("deviceid");
        const std::size_t i_category                     = get_index("category");
        const std::size_t i_software_first               = get_index("softwarefirst");
        const std::size_t i_embodied_co2_kg              = get_index("embodiedco2_kg");
        const std::size_t i_expected_lifetime_years      = get_index("expectedlifetime_years");
        const std::size_t i_annual_energy_savings_kwh    = get_index("annual_energy_savings_kwh");
        const std::size_t i_annual_pollutant_reduction   = get_index("annual_pollutant_reduction_tons");
        const std::size_t i_ecoimpact_score              = get_index("ecoimpactscore_01");
        const std::size_t i_notes                        = get_index("notes");

        std::vector<EcoDeviceRecord> records;
        std::string line;
        while (std::getline(file, line))
        {
            if (line.empty())
            {
                continue;
            }

            const auto cols = detail::split_csv_line(line);
            if (cols.size() < header_cols.size())
            {
                continue;
            }

            EcoDeviceRecord rec;
            rec.device_id                     = cols[i_device_id];
            rec.category                      = cols[i_category];
            rec.software_first                = detail::parse_bool(cols[i_software_first]);
            rec.embodied_co2_kg               = detail::parse_double(cols[i_embodied_co2_kg]);
            rec.expected_lifetime_years       = detail::parse_double(cols[i_expected_lifetime_years]);
            rec.annual_energy_savings_kwh     = detail::parse_double(cols[i_annual_energy_savings_kwh]);
            rec.annual_pollutant_reduction_tons = detail::parse_double(cols[i_annual_pollutant_reduction_tons]);
            rec.ecoimpact_score_01            = detail::parse_double(cols[i_ecoimpact_score]);
            rec.notes                         = cols[i_notes];

            records.push_back(rec);
        }

        return records;
    }

    inline EcoDeviceCorridorAggregate aggregate_devices_simple(
        const std::vector<EcoDeviceRecord>& devices,
        const std::string& corridor_id = "GLOBAL"
    )
    {
        EcoDeviceCorridorAggregate agg;
        agg.corridor_id = corridor_id;

        double ecoimpact_sum = 0.0;

        for (const auto& d : devices)
        {
            agg.total_embodied_co2_kg               += d.embodied_co2_kg;
            agg.total_annual_energy_savings_kwh     += d.annual_energy_savings_kwh;
            agg.total_annual_pollutant_reduction_tons += d.annual_pollutant_reduction_tons;
            ecoimpact_sum                           += d.ecoimpact_score_01;
            agg.device_count                        += 1;
        }

        if (agg.device_count > 0)
        {
            agg.average_ecoimpact_score_01 = ecoimpact_sum / static_cast<double>(agg.device_count);
        }
        else
        {
            agg.average_ecoimpact_score_01 = 0.0;
        }

        return agg;
    }

    inline double compute_device_lifecycle_score(
        const EcoDeviceRecord& rec,
        double max_embodied_co2_kg,
        double max_energy_savings_kwh,
        double max_pollutant_reduction_tons
    )
    {
        const double alpha = 0.4;
        const double beta  = 0.3;
        const double gamma = 0.3;

        double embodied_term = 0.0;
        if (max_embodied_co2_kg > 0.0)
        {
            const double norm_embodied = std::min(rec.embodied_co2_kg / max_embodied_co2_kg, 1.0);
            embodied_term = 1.0 - norm_embodied;
        }

        double energy_term = 0.0;
        if (max_energy_savings_kwh > 0.0)
        {
            const double norm_energy = std::min(rec.annual_energy_savings_kwh / max_energy_savings_kwh, 1.0);
            energy_term = norm_energy;
        }

        double pollutant_term = 0.0;
        if (max_pollutant_reduction_tons > 0.0)
        {
            const double norm_poll = std::min(rec.annual_pollutant_reduction_tons / max_pollutant_reduction_tons, 1.0);
            pollutant_term = norm_poll;
        }

        double score = alpha * embodied_term +
                       beta  * energy_term +
                       gamma * pollutant_term;

        score = std::max(0.0, std::min(score, 1.0));
        return score;
    }

    inline void compute_lifecycle_scores_for_all(
        std::vector<EcoDeviceRecord>& devices,
        std::vector<double>& out_scores
    )
    {
        out_scores.clear();
        out_scores.reserve(devices.size());

        double max_embodied = 0.0;
        double max_energy   = 0.0;
        double max_poll     = 0.0;

        for (const auto& d : devices)
        {
            if (d.embodied_co2_kg > max_embodied) max_embodied = d.embodied_co2_kg;
            if (d.annual_energy_savings_kwh > max_energy) max_energy = d.annual_energy_savings_kwh;
            if (d.annual_pollutant_reduction_tons > max_poll) max_poll = d.annual_pollutant_reduction_tons;
        }

        for (const auto& d : devices)
        {
            const double s = compute_device_lifecycle_score(
                d,
                max_embodied,
                max_energy,
                max_poll
            );
            out_scores.push_back(s);
        }
    }

} // namespace ceim

#endif // CEIM_ECODEVICESCORE_HPP
