"""Implementation of grain-assisted recombination"""

# MyFloat alpha_recomb_grain(int i, MyFloat temp, MyFloat x_elec, MyFloat shieldfac, char *ion_name)
# {
#     MyFloat psi = grain_charge_psi(i, temp, x_elec, shieldfac);
#     // MyFloat temp = get_temperature(i);
#     int j = ion_name_to_index(ion_name);

#     MyFloat C[NUM_RECOMB_TABLE_IONS][7] = {
#         {12.25, 8.074E-6, 1.378, 5.087E2, 1.586E-2, 0.4723, 1.102E-5}, // H+
#         {5.572, 3.185E-7, 1.512, 5.115E3, 3.903E-7, 0.4956, 5.494E-7}, // He+
#         {45.58, 6.089E-3, 1.128, 4.331E2, 4.845E-2, 0.8120, 1.333E-4}, // C+
#         {2.178, 1.732E-7, 2.133, 1.029E4, 1.859E-6, 1.0341, 3.223E-5}, // Na+
#         {2.510, 8.116E-8, 1.864, 6.170E4, 2.169E-6, 0.9605, 7.232E-5}, // Mg+
#         {2.166, 5.678E-8, 1.874, 4.375E4, 1.635E-6, 0.8964, 7.538E-5}, // Si+
#         {3.064, 7.769E-5, 1.319, 1.087E2, 3.475E-1, 0.4790, 4.689E-2}, // S+
#         {1.596, 1.907E-7, 2.123, 8.138E3, 1.530E-5, 1.0380, 4.550E-5}, // K+
#         {1.636, 8.208E-9, 2.289, 1.254E5, 1.349E-9, 1.1506, 7.204E-4}, // Ca+
#         {2.029, 1.433E-6, 1.673, 1.403E4, 1.865E-6, 0.9358, 4.339E-9}, // Mn+
#         {1.701, 9.554E-8, 1.851, 5.763E4, 4.116E-8, 0.9456, 2.198E-5}, // Fe+
#         {8.270, 2.051E-4, 1.252, 1.590E2, 6.072E-2, 0.5980, 4.497E-7}  // Ca++
#     };

#     MyFloat Z = P[i].Metallicity[0] / All.SolarAbundances[0];
#     return Z * 1e-14 * C[j][0] / (1 + C[j][1] * pow(psi, C[j][2]) * (1 + C[j][3] * pow(temp, C[j][4]) * pow(psi, -C[j][5] - C[j][6] * log(temp))));
# }
