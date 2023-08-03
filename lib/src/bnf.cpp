#include "bnf.h"

void bnf(bool &disp, unsigned short int &n, unsigned int &s, vector<double> &excess, double &prec, vector<bool> &unsettled, unsigned short int &iter, unsigned int &piv, unsigned int &sr, double &t, vector<double> &x, vector<vector<bool>> &A, double &t1, vector<double> &singleton_bounds, bool &nlsu)
{
	vector<bool> unsettled_p(n, true);
	vector<vector<double>> Arref(n, vector<double>(n, 0));
	Arref[0] = vector<double>(n, 1);
	vector<bool> J(n, true);
	J[0] = false;
	unsigned short int rank = 1;
	vector<vector<bool>> Asettled(n, vector<bool>(n, 0));
	Asettled[0] = vector<bool>(n, true);
	if (disp)
	{
		cout << endl
			 << "   ---===   ITERATION " << iter + 1 << "   ===---   " << endl
			 << endl;
		cout << "Starting point:" << endl;
		for (unsigned short int i = 0; i < n; i++)
		{
			cout << x[i] << endl;
		}
	}
	
	vector<double> d(n, 0);
	double epsi = 0;
	double epsi_old = -DBL_MAX;
	while (rank < n)
	{
		pivot(epsi, s, excess, prec, n, A, Arref, J, unsettled, rank, d, x, disp, Asettled, piv, sr, iter, unsettled_p, singleton_bounds, epsi_old, nlsu);
	}
	t = cpuTime() - t1;
	if (disp)
	{
		cout << "BNF finished!" << endl;
		cout << "The nucleolus solution:" << endl;

		for (unsigned short int i = 0; i < n; i++)
		{
			cout << x[i] << endl;
		}

		cout << "Time needed: " << t << " seconds" << endl;
		cout << "Iterations needed: " << iter << endl;
		cout << "Pivots needed: " << piv << endl;
		cout << "Subroutine solves needed: " << sr << endl;
	}
}

void pivot(double &epsi, unsigned int &s, vector<double> &excess, double &prec, unsigned short int &n, vector<vector<bool>> &A, vector<vector<double>> &Arref, vector<bool> &J, vector<bool> &unsettled, unsigned short int &rank, vector<double> &d, vector<double> &x, bool &disp, vector<vector<bool>> &Asettled, unsigned int &piv, unsigned int &sr_count, unsigned short int &iter, vector<bool> &unsettled_p, vector<double> &singleton_bounds, double &epsi_old, bool &nlsu)
{
	vec_min_uns(epsi, excess, unsettled, s);
	if (disp)
	{
		cout << "Epsilon: " << epsi << endl;
	}
	vector<bool> T(s, false);
	vector<unsigned int> T_coord(0, 0);
	vector<bool> T2(n, false);
	vector<unsigned int> T2_coord(0, 0);
	unsigned int t_size = 0;
	tight_coal(T, excess, epsi, prec, s, T_coord, unsettled, t_size);
	unsigned short int t2_size = 0;
	tight_coal2(T2, x, singleton_bounds, prec, n, T2_coord, unsettled_p, t2_size);
	vector<vector<bool>> Atight(t_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t_size; i++)
	{
		Atight[i] = A[T_coord[i]];
	}
	vector<vector<bool>> Atight2(t2_size, vector<bool>(n, false));
	for (unsigned int i = 0; i < t2_size; i++)
	{
		Atight2[i] = A[T2_coord[i]];
	}
	vector<bool> U(t_size, true);
	vector<bool> U2(t2_size, true);
	bool u = true;
	bool settled = false;
	subroutine(U, U2, Atight, Atight2, Arref, J, prec, n, t_size, t2_size, rank, disp, Asettled, sr_count, u, s, T_coord, T2_coord, unsettled, epsi_old, epsi, unsettled_p, settled, nlsu);
	if (disp)
	{
		cout << endl
			 << "  --==  subroutine finished  ==--  " << endl
			 << endl;
	}
	if (settled)
	{
		iter++;
	}
	if (disp)
	{
		cout << "T:" << endl;
		for (unsigned int i = 0; i < t_size; i++)
		{
			if (!U[i])
			{
				cout << T_coord[i] + 1 << endl;
			}
		}
		cout << "U:" << endl;
		for (unsigned int i = 0; i < t_size; i++)
		{
			if (U[i])
			{
				cout << T_coord[i] + 1 << endl;
			}
		}
		cout << "T0:" << endl;
		for (unsigned int i = 0; i < t2_size; i++)
		{
			if (!U2[i])
			{
				cout << T2_coord[i] + 1 << endl;
			}
		}
		cout << "U0:" << endl;
		for (unsigned int i = 0; i < t2_size; i++)
		{
			if (U2[i])
			{
				cout << T2_coord[i] + 1 << endl;
			}
		}
	}
	if (u)
	{
		piv++;
		if (disp)
		{
			cout << endl
				 << "  --==  solving improving direction LP  ==--  " << endl
				 << endl;
		}
		imprdir(d, n, t_size, t2_size, Atight, Atight2, U, U2, rank, Asettled, disp);
		if (disp)
		{
			cout << endl
				 << "  --==  improving direction obtained  ==--  " << endl
				 << endl;
		}
		if (disp)
		{
			cout << "Improving direction:" << endl;
			for (unsigned short int i = 0; i < n; i++)
			{
				cout << d[i] << "    ";
			}
			cout << endl;
		}
		if (disp)
		{
			cout << endl
				 << "  --==  computing step size  ==--  " << endl
				 << endl;
		}
		step(T, T2, unsettled, unsettled_p, s, A, epsi, excess, d, n, x, singleton_bounds, disp, prec);
	}
	else
	{
		if (disp)
		{
			cout << endl
				 << "  --==  minimal tight set found!  ==--  " << endl
				 << endl;
		}
		if (rank == n)
		{
			return;
		}
		if (!nlsu)
		{
			for (unsigned int i = 0; i < s; i++)
			{
				if (unsettled[i])
				{
					if (!(binrank(Arref, J, A[i], n)))
					{
						unsettled[i] = false;
						unsettled[s - 1 - i] = false;
					}
				}
			}
		}
		for (unsigned short int i = 0; i < n; i++)
		{
			if (unsettled_p[i] == true && unsettled[pow(2, i) - 1] == false)
			{
				unsettled_p[i] = false;
			}
		}
	}
	if (disp && settled)
	{
		cout << endl
			 << "   ---===   ITERATION " << iter + 1 << "   ===---   " << endl
			 << endl;
	}
}

void subroutine(vector<bool> &U, vector<bool> &U2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<vector<double>> &Arref, vector<bool> &J, double &prec, unsigned short int &n, unsigned int &tight_size, unsigned short int &tight2_size, unsigned short int &rank, bool &disp, vector<vector<bool>> &Asettled, unsigned int &sr_count, bool &u, unsigned int &s, vector<unsigned int> &T_coord, vector<unsigned int> &T2_coord, vector<bool> &unsettled, double &epsi_old, double &epsi, vector<bool> &unsettled_p, bool &settled, bool &nlsu)
{
	unsigned int sumt = 0;
	vector<bool> t(tight_size, false);
	unsigned int sumt2 = 0;
	vector<bool> t2(tight2_size, false);

	glp_term_out(GLP_OFF);
	glp_prob *lp;
	int ia[1 + 1000] = {0};
	int ja[1 + 1000] = {0};
	double ar[1 + 1000] = {0};
	lp = glp_create_prob();
	glp_set_prob_name(lp, "P(1)");
	glp_set_obj_dir(lp, GLP_MAX);
	glp_add_cols(lp, tight_size + tight2_size + rank);

	for (unsigned int i = 0; i < tight_size + tight2_size + rank; i++)
	{
		const int col_index = i + 1;
		const std::string col_name_string = "lambda{" + std::to_string(col_index) + "}";
		const char *col_name = col_name_string.c_str();
		glp_set_col_name(lp, col_index, col_name);

		if (i < tight_size + tight2_size)
		{
			glp_set_col_bnds(lp, col_index, GLP_LO, 0.0, 0.0);
			glp_set_obj_coef(lp, col_index, 1);
		}
		else
		{
			glp_set_col_bnds(lp, col_index, GLP_FR, 0.0, 0.0);
			glp_set_obj_coef(lp, col_index, 0.0);
		}
	}

	int ia_index = 1;
	vector<int> q_indices(tight_size + tight2_size + rank, -1);
	int constr_index = 0;

	constr_index++;
	glp_add_rows(lp, 1);
	glp_set_row_name(lp, constr_index, "q = 1");
	glp_set_row_bnds(lp, constr_index, GLP_FX, 1, 0.0);
	int q_constr_index = constr_index;

	for (unsigned int j = 0; j < tight_size; j++)
	{
		const int col_index = j + 1;
		q_indices[j] = ia_index;
		ia[ia_index] = constr_index;
		ja[ia_index] = col_index;
		ar[ia_index] = 1.0;
		ia_index++;
	}

	for (unsigned short int i = 0; i < n; i++)
	{
		constr_index++;
		glp_add_rows(lp, 1);
		const std::string row_name_string = "p(" + std::to_string(i) + ") = 0";
		const char *row_name = row_name_string.c_str();
		glp_set_row_name(lp, constr_index, row_name);
		glp_set_row_bnds(lp, constr_index, GLP_FX, 0.0, 0.0);

		for (unsigned int j = 0; j < tight_size; j++)
		{
			if (Atight[j][i])
			{
				const int col_index = j + 1;
				ia[ia_index] = constr_index;
				ja[ia_index] = col_index;
				ar[ia_index] = 1.0;
				ia_index++;
			}
			if (j < rank && Asettled[j][i])
			{
				const int col_index = j + tight_size + tight2_size + 1;
				ia[ia_index] = constr_index;
				ja[ia_index] = col_index;
				ar[ia_index] = 1.0;
				ia_index++;
			}
		}
		for (unsigned int j = 0; j < tight2_size; j++)
		{
			if (Atight2[j][i])
			{
				const int col_index = j + tight_size + 1;
				ia[ia_index] = constr_index;
				ja[ia_index] = col_index;
				ar[ia_index] = 1.0;
				ia_index++;
			}
		}
		if (rank > tight_size)
		{
			for (unsigned short int j = tight_size; j < rank; j++)
			{
				if (Asettled[j][i])
				{
					const int col_index = j + tight_size + tight2_size + 1;
					ia[ia_index] = constr_index;
					ja[ia_index] = col_index;
					ar[ia_index] = 1.0;
					ia_index++;
				}
			}
		}
	}

	glp_load_matrix(lp, ia_index - 1, ia, ja, ar);
	glp_simplex(lp, nullptr);

	int lp_status = glp_get_status(lp);
	bool lp_optimal = lp_status == GLP_OPT;
	bool lp_feas = lp_status == GLP_FEAS;
	bool feas = lp_optimal || lp_feas;

	if (disp)
	{
		cout << endl
			 << "  --==  solving subroutine LP  ==--  " << endl
			 << endl;
	}

	if (disp)
	{
		cout << "subroutine feasibility: " << feas << endl;
	}
	if (feas && nlsu)
	{
		settled = true;
	}

	sr_count++;
	unsigned int i;
	while (feas)
	{
		subr_upd(Arref, J, i, n, prec, U, U2, sumt, sumt2, t, t2, Atight, Atight2, tight_size, tight2_size, rank, unsettled, Asettled, disp, s, T_coord, T2_coord, epsi_old, epsi, unsettled_p, settled, lp, ia, ja, ar, ia_index, q_constr_index, q_indices);
		if (rank == n)
		{
			u = false;
			return;
		}
		if (sumt < tight_size)
		{
			i = 0;
			while (i < tight_size)
			{
				if (!t[i])
				{
					if (!(binrank(Arref, J, Atight[i], n)))
					{
						U[i] = false;
						t[i] = true;

						int ia_ref_index = q_indices[i];
						if (ia_ref_index > -1)
						{
							ar[ia_ref_index] = ar[ia_ref_index] - 1;
						}
						else
						{
							const int col_index = i + 1;
							q_indices[i] = ia_index;
							ia[ia_index] = q_constr_index;
							ja[ia_index] = col_index;
							ar[ia_index] = -1;
							ia_index++;
						}

						const int obj_index = i + 1;
						const double sr_obj_coef = glp_get_obj_coef(lp, obj_index);
						glp_set_obj_coef(lp, obj_index, sr_obj_coef - 1);

						sumt++;
						unsettled[T_coord[i]] = false;
						unsettled[s - 1 - T_coord[i]] = false;
						if (disp)
						{
							cout << T_coord[i] + 1 << " and " << s - T_coord[i] << " got settled without rank increase."
								 << endl;
						}
						if (sumt == tight_size && sumt2 == tight2_size)
						{
							u = false;
							return;
						}
					}
				}
				i++;
			}
			i = 0;
			while (i < tight2_size)
			{
				if (!t2[i])
				{
					if (!(binrank(Arref, J, Atight2[i], n)))
					{
						U2[i] = false;
						t2[i] = true;

						const int obj_index = i + tight_size + 1;
						const double sr_obj_coef = glp_get_obj_coef(lp, obj_index);
						glp_set_obj_coef(lp, obj_index, sr_obj_coef - 1);

						sumt2++;
						unsettled[T2_coord[i]] = false;
						unsettled[s - 1 - T2_coord[i]] = false;
						if (disp)
						{
							cout << T2_coord[i] + 1 << " and " << s - T2_coord[i]
								 << " got settled without rank increase." << endl;
						}
						if (sumt == tight_size && sumt2 == tight2_size)
						{
							u = false;
							return;
						}
					}
				}
				i++;
			}
			for (unsigned short int i = 0; i < n; i++)
			{
				if (unsettled_p[i] && !unsettled[pow(2, i) - 1])
				{
					unsettled_p[i] = false;
				}
			}

			glp_set_obj_dir(lp, GLP_MIN);
			glp_load_matrix(lp, ia_index - 1, ia, ja, ar);
			glp_simplex(lp, nullptr);

			lp_status = glp_get_status(lp);
			lp_optimal = lp_status == GLP_OPT;
			lp_feas = lp_status == GLP_FEAS;
			feas = lp_optimal || lp_feas;

			if (disp)
			{
				cout << endl
					 << "  --==  solving subroutine LP again  ==--  " << endl
					 << endl;
			}

			if (disp)
			{
				cout << "subroutine feasibility: " << feas << endl;
			}
			sr_count++;
		}
		else
		{
			u = false;
			return;
		}
	}

	glp_delete_prob(lp);
	return;
}

void subr_upd(vector<vector<double>> &Arref, vector<bool> &J, unsigned int &i, unsigned short int &n, double &prec, vector<bool> &U, vector<bool> &U2, unsigned int &sumt, unsigned int &sumt2, vector<bool> &t, vector<bool> &t2, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, unsigned int &tight_size, unsigned short int &tight2_size, unsigned short int &rank, vector<bool> &unsettled, vector<vector<bool>> &Asettled, bool &disp, unsigned int &s, vector<unsigned int> &T_coord, vector<unsigned int> &T2_coord, double &epsi_old, double &epsi, vector<bool> &unsettled_p, bool &settled, glp_prob *lp, int ia[], int ja[], double ar[], int &ia_index, int &q_constr_index, vector<int> &q_indices)
{
	i = 0;
	vector<double> lambdi(tight_size + tight2_size, 0);
	for (unsigned int j = 0; j < tight_size; j++)
	{
		if (!t[j])
		{
			const int col_index = j + 1;
			lambdi[j] = glp_get_col_prim(lp, col_index);
		}
	}
	for (unsigned short int j = 0; j < tight2_size; j++)
	{
		if (!t2[j])
		{
			const int col_index = j + tight_size + 1;
			lambdi[j + tight_size] = glp_get_col_prim(lp, col_index);
		}
	}

	while (i < tight_size && sumt < tight_size)
	{
		if (lambdi[i] > prec)
		{
			U[i] = false;
			t[i] = true;

			int ia_ref_index = q_indices[i];
			if (ia_ref_index > -1)
			{
				ar[ia_ref_index] = ar[ia_ref_index] - 1;
			}
			else
			{
				const int col_index = i + 1;
				q_indices[i] = ia_index;
				ia[ia_index] = q_constr_index;
				ja[ia_index] = col_index;
				ar[ia_index] = -1;
				ia_index++;
			}

			const int obj_index = i + 1;
			const double sr_obj_coef = glp_get_obj_coef(lp, obj_index);
			glp_set_obj_coef(lp, obj_index, sr_obj_coef - 1);

			sumt++;
			unsettled[T_coord[i]] = false;
			unsettled[s - 1 - T_coord[i]] = false;
			if (binrank(Arref, J, Atight[i], n))
			{
				rank++;
				if (epsi > epsi_old)
				{
					settled = true;
					epsi_old = epsi;
				}
				if (disp)
				{
					cout << "lambda_" << T_coord[i] + 1 << " > 0, rank = " << rank << " (" << s - T_coord[i]
						 << " settled as well)" << endl;
				}
				rowechform(Arref, J, Atight[i], n, rank);
				Asettled[rank - 1] = Atight[i];
				if (rank == n)
				{
					if (disp)
					{
						cout << "Rank condition satisfied!" << endl;
					}
					return;
				}

				const int col_index = i + 1;
				glp_set_col_bnds(lp, col_index, GLP_FR, 0.0, 0.0);
			}
			else
			{
				if (disp)
				{
					cout << "lambda_" << T_coord[i] + 1 << " > 0, got settled (with " << s - T_coord[i]
						 << ") without rank increase" << endl;
				}
			}
		}
		i++;
	}
	i = 0;
	while (i < tight2_size && sumt2 < tight2_size)
	{
		if (lambdi[i + tight_size] > prec)
		{
			U2[i] = false;
			t2[i] = true;
			sumt2++;

			const int obj_index = i + tight_size + 1;
			const double sr_obj_coef = glp_get_obj_coef(lp, obj_index);
			glp_set_obj_coef(lp, obj_index, sr_obj_coef - 1);

			unsettled[T2_coord[i]] = false;
			unsettled[s - 1 - T2_coord[i]] = false;
			if (binrank(Arref, J, Atight2[i], n))
			{
				rank++;
				if (epsi > epsi_old)
				{
					settled = true;
					epsi_old = epsi;
				}
				if (disp)
				{
					cout << "lambda_" << T2_coord[i] + 1 << " > 0, rank = " << rank << " (" << s - T2_coord[i]
						 << " settled as well)" << endl;
				}
				rowechform(Arref, J, Atight2[i], n, rank);
				Asettled[rank - 1] = Atight2[i];
				if (rank == n)
				{
					if (disp)
					{
						cout << "Rank condition satisfied!" << endl;
					}
					return;
				}

				const int col_index = i + tight_size + 1;
				glp_set_col_bnds(lp, col_index, GLP_FR, 0.0, 0.0);
			}
			else
			{
				if (disp)
				{
					cout << "lambda_" << T_coord[i] + 1 << " > 0, got settled (with " << s - T_coord[i]
						 << ") without rank increase" << endl;
				}
			}
		}
		i++;
	}
	for (unsigned short int i = 0; i < n; i++)
	{
		if (unsettled_p[i] && !unsettled[pow(2, i) - 1])
		{
			unsettled_p[i] = false;
		}
	}
}

void imprdir(vector<double> &d, unsigned short int &n, unsigned int &t_size, unsigned short int &t2_size, vector<vector<bool>> &Atight, vector<vector<bool>> &Atight2, vector<bool> &U, vector<bool> &U2, unsigned short int &rank, vector<vector<bool>> &Asettled, bool &disp)
{
	glp_term_out(GLP_OFF);
	glp_prob *lp;
	int ia[1 + 1000] = {0};
	int ja[1 + 1000] = {0};
	double ar[1 + 1000] = {0};
	lp = glp_create_prob();
	glp_set_prob_name(lp, "IMPRDIR");
	glp_set_obj_dir(lp, GLP_MIN);
	glp_add_cols(lp, n);

	for (unsigned int i = 0; i < n; i++)
	{
		const int col_index = i + 1;
		const std::string col_name_string = "dir(" + std::to_string(col_index) + ")";
		const char *col_name = col_name_string.c_str();
		glp_set_col_name(lp, col_index, col_name);
		glp_set_col_bnds(lp, col_index, GLP_FR, 0.0, 0.0);
	}

	int ia_index = 1;
	int constr_index = 0;
	for (unsigned int i = 0; i < t_size; i++)
	{
		constr_index++;
		glp_add_rows(lp, 1);
		const std::string row_name_string = "ineq_t_size(" + std::to_string(i) + ")";
		const char *row_name = row_name_string.c_str();
		glp_set_row_name(lp, constr_index, row_name);

		for (unsigned short int j = 0; j < n; j++)
		{
			if (Atight[i][j])
			{
				const int col_index = j + 1;
				glp_set_obj_coef(lp, col_index, 1);
				ia[ia_index] = constr_index;
				ja[ia_index] = col_index;
				ar[ia_index] = 1.0;
				ia_index++;
			}
		}

		if (U[i])
		{
			glp_set_row_bnds(lp, constr_index, GLP_LO, 1, 0.0);
		}
		else
		{
			glp_set_row_bnds(lp, constr_index, GLP_FX, 0.0, 0.0);
		}
	}

	for (unsigned int i = 0; i < t2_size; i++)
	{
		constr_index++;
		glp_add_rows(lp, 1);
		const std::string row_name_string = "ineq_t2_size(" + std::to_string(i) + ") >= 0";
		const char *row_name = row_name_string.c_str();
		glp_set_row_name(lp, constr_index, row_name);
		glp_set_row_bnds(lp, constr_index, GLP_LO, 0.0, 0.0);

		for (unsigned short int j = 0; j < n; j++)
		{
			if (Atight2[i][j])
			{
				const int col_index = j + 1;
				ia[ia_index] = constr_index;
				ja[ia_index] = col_index;
				ar[ia_index] = 1.0;
				ia_index++;
			}
		}
	}

	for (unsigned short int i = 0; i < rank; i++)
	{
		constr_index++;
		glp_add_rows(lp, 1);
		const std::string row_name_string = "eq_rank(" + std::to_string(i) + ") = 0";
		const char *row_name = row_name_string.c_str();
		glp_set_row_name(lp, constr_index, row_name);
		glp_set_row_bnds(lp, constr_index, GLP_FX, 0.0, 0.0);

		for (unsigned short int j = 0; j < n; j++)
		{
			if (Asettled[i][j])
			{
				const int col_index = j + 1;
				ia[ia_index] = constr_index;
				ja[ia_index] = col_index;
				ar[ia_index] = 1.0;
				ia_index++;
			}
		}
	}

	glp_load_matrix(lp, ia_index - 1, ia, ja, ar);
	glp_simplex(lp, nullptr);

	for (unsigned short int i = 0; i < n; i++)
	{
		const int col_index = i + 1;
		d[i] = glp_get_col_prim(lp, col_index);
	}

	glp_delete_prob(lp);
	return;
}

void step(vector<bool> &T, vector<bool> &T2, vector<bool> &unsettled, vector<bool> &unsettled_p, unsigned int &s, vector<vector<bool>> &A, double &epsi, vector<double> &excess, vector<double> &d, unsigned short int &n, vector<double> &x, vector<double> &singleton_bounds, bool &disp, double &prec)
{
	double alpha = DBL_MAX;
	double Ad;
	for (unsigned int i = 0; i < s; i++)
	{
		if (!T[i] && unsettled[i])
		{
			Ad = 0;
			for (unsigned short int j = 0; j < n; j++)
			{
				if (A[i][j])
					Ad += d[j];
			}
			if (Ad < 1 - prec && (epsi - excess[i]) / (Ad - 1) < alpha)
			{
				alpha = (epsi - excess[i]) / (Ad - 1);
			}
		}
	}
	for (unsigned short int i = 0; i < n; i++)
	{
		if (!T2[i] && unsettled_p[i])
		{
			if (d[i] < -prec && (singleton_bounds[i] - x[i]) / d[i] < alpha)
			{
				alpha = (singleton_bounds[i] - x[i]) / d[i];
			}
		}
	}
	if (disp)
	{
		cout << "Step size: " << alpha << endl;
	}
	if (disp)
	{
		cout << endl
			 << "  --==  step size obtained  ==--  " << endl
			 << endl;
	}
	for (unsigned short int i = 0; i < n; i++)
	{
		x[i] += alpha * d[i];
	}
	if (disp)
	{
		cout << "New point: " << endl;
		for (unsigned short int i = 0; i < n; i++)
		{
			cout << x[i] << endl;
		}
	}
	for (unsigned int i = 0; i < s; i++)
	{
		if (unsettled[i])
		{
			Ad = 0;
			for (unsigned short int j = 0; j < n; j++)
			{
				if (A[i][j])
				{
					Ad += d[j];
				}
			}
			excess[i] += alpha * Ad;
		}
	}
}

void tight_coal(vector<bool> &T, vector<double> &excess, double &epsi, double &prec, unsigned int &s, vector<unsigned int> &T_coord, vector<bool> &unsettled, unsigned int &t_size)
{
	for (unsigned int i = 0; i < s; i++)
	{
		if (unsettled[i])
		{
			if (abs(excess[i] - epsi) < prec)
			{
				t_size++;
				T[i] = true;
				T_coord.push_back(i);
			}
		}
	}
}

void tight_coal2(vector<bool> &T2, vector<double> &x, vector<double> &singleton_bounds, double &prec, unsigned short int &n, vector<unsigned int> &T2_coord, vector<bool> &unsettled_p, unsigned short int &t2_size)
{
	for (unsigned int i = 0; i < n; i++)
	{
		if (unsettled_p[i])
		{
			if (abs(x[i] - singleton_bounds[i]) < prec)
			{
				t2_size++;
				T2[i] = true;
				T2_coord.push_back(pow(2, i) - 1);
			}
		}
	}
}
