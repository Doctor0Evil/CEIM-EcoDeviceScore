#include <iostream>
#include <iomanip>
#include <string>
#include "ceim_ecodevice_score.hpp"

int main()
{
    try
    {
        const std::string csv_path =
            "../qpudatashards/particles/EcoTechnologyDeviceLifecycleScore2026v1.csv";

        auto devices = ceim::load_eco_device_csv(csv_path);

        auto agg = ceim::aggregate_devices_simple(devices, "PHX-CORRIDOR-ALL");

        std::cout << "Corridor ID: " << agg.corridor_id << "\n";
        std::cout << "Device count: " << agg.device_count << "\n";
        std::cout << "Total embodied CO2 [kg]: " << agg.total_embodied_co2_kg << "\n";
        std::cout << "Total annual energy savings [kWh]: "
                  << agg.total_annual_energy_savings_kwh << "\n";
        std::cout << "Total annual pollutant reduction [tons]: "
                  << agg.total_annual_pollutant_reduction_tons << "\n";
        std::cout << "Average ecoimpact score [0-1]: "
                  << std::fixed << std::setprecision(3)
                  << agg.average_ecoimpact_score_01 << "\n\n";

        std::vector<double> lifecycle_scores;
        ceim::compute_lifecycle_scores_for_all(devices, lifecycle_scores);

        std::cout << "Per-device lifecycle scores:\n";
        for (std::size_t i = 0; i < devices.size(); ++i)
        {
            const auto& d = devices[i];
            const double s = lifecycle_scores[i];

            std::cout << "  " << d.device_id
                      << " (" << d.category << ")"
                      << " -> lifecycle_score=" << std::fixed << std::setprecision(3)
                      << s
                      << ", ecoimpact_score_01=" << d.ecoimpact_score_01
                      << "\n";
        }

        std::vector<ceim::CeimNodeSample> samples;
        samples.push_back(ceim::CeimNodeSample{
            5.0,   // Cin
            2.0,   // Cout
            1.0,   // Q (m^3/s)
            86400, // dt (1 day)
            10.0,  // Cref
            1.0    // lambda
        });

        const auto kn_result = ceim::compute_node_impact_with_uncertainty(samples);
        std::cout << "\nExample CEIM node impact (single sample):\n";
        std::cout << "  Kn = " << kn_result.kn
                  << " (std ~ " << kn_result.kn_std << ")\n";
    }
    catch (const std::exception& ex)
    {
        std::cerr << "Error: " << ex.what() << "\n";
        return 1;
    }

    return 0;
}
