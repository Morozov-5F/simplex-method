#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define PRINT_FORMAT_ERROR            "\x1B[31m"
#define PRINT_FORMAT_HIGHLIGHT        "\x1B[35m"
#define PRINT_FORMAT_HIGHLIGHT_2      "\x1B[32m"
#define PRINT_CLEAR_FORMAT            "\x1B[0m"

#define FILE_NAME                     "input.txt"

void dump_matrix(gsl_matrix * M)
{
   size_t m = M->size1, n = M->size2;

   for (int i = 0; i < m; ++i)
   {
      for (int j = 0; j < n; ++j)
      {
         double elem = gsl_matrix_get(M, i, j);
         printf("% 03.2f\t", elem);
      }
      printf("\n");
   }
   printf("\n");
}

void dump_matrix_h_column(gsl_matrix * M, size_t column)
{
   size_t m = M->size1, n = M->size2;

   for (int i = 0; i < m; ++i)
   {
      for (int j = 0; j < n; ++j)
      {
         double elem = gsl_matrix_get(M, i, j);
         if (j == column)
         {
            printf(PRINT_FORMAT_HIGHLIGHT "% 03.2f\t" PRINT_CLEAR_FORMAT, elem);
         }
         else
         {
            printf("% 03.2f\t", elem);
         }
      }
      printf("\n");
   }
   printf("\n");
}

void dump_matrix_h_row(gsl_matrix * M, size_t row)
{
   size_t m = M->size1, n = M->size2;

   for (int i = 0; i < m; ++i)
   {
      for (int j = 0; j < n; ++j)
      {
         double elem = gsl_matrix_get(M, i, j);
         if (i == row)
         {
            printf(PRINT_FORMAT_HIGHLIGHT "% 03.2f\t" PRINT_CLEAR_FORMAT, elem);
         }
         else
         {
            printf("% 03.2f\t", elem);
         }
      }
      printf("\n");
   }
   printf("\n");
}

void dump_matrix_h_elem(gsl_matrix * M, size_t column, size_t row)
{
   size_t m = M->size1, n = M->size2;

   for (int i = 0; i < m; ++i)
   {
      for (int j = 0; j < n; ++j)
      {
         double elem = gsl_matrix_get(M, i, j);
         if (j == column && i == row)
         {
            printf(PRINT_FORMAT_HIGHLIGHT "% 03.2f\t" PRINT_CLEAR_FORMAT, elem);
         }
         else
         {
            printf("% 03.2f\t", elem);
         }
      }
      printf("\n");
   }
}

void dump_matrix_result(gsl_matrix * M)
{
   size_t m = M->size1, n = M->size2;

   for (int i = 0; i < m; ++i)
   {
      for (int j = 0; j < n; ++j)
      {
         double elem = gsl_matrix_get(M, i, j);
         if (j == n - 1)
         {
            printf(PRINT_FORMAT_HIGHLIGHT_2 "% 03.2f\t" PRINT_CLEAR_FORMAT, elem);
         }
         else
         {
            printf("% 03.2f\t", elem);
         }
      }
      printf("\n");
   }
}

int main(int argc, const char * argv[])
{
   unsigned m, n;
   size_t lead_index_column, lead_index_row;
   double lead_elem;
   gsl_matrix * matrix;
   gsl_vector * lead_column, * res_column, * lead_row_copy;
   gsl_vector_view z_row, lead_row, current_row;
   gsl_matrix_view main_matrix;

   FILE * input_file = fopen(FILE_NAME, "r");

   if (!input_file)
   {
      fprintf(stderr, PRINT_FORMAT_ERROR "Error has occured while openning file %s.\n" PRINT_CLEAR_FORMAT, FILE_NAME);
      goto done;
   }

   fscanf(input_file, "%u %u", &m, &n);
   lead_column = gsl_vector_alloc(m);
   res_column = gsl_vector_alloc(m);
   lead_row_copy = gsl_vector_alloc(n);
   matrix = gsl_matrix_alloc(m, n);
   gsl_matrix_fscanf(input_file, matrix);

   fclose(input_file);

   /* Simplex algorithm */
   z_row = gsl_matrix_subrow(matrix, m - 1, 0, m - 1);
   main_matrix = gsl_matrix_submatrix(matrix, 0, 0, m - 1, n);

   dump_matrix(matrix);

   /* While we have negative elements in last row */
   while (gsl_vector_get(&z_row.vector,
          lead_index_column = gsl_vector_min_index(&z_row.vector)) < 0)
   {
      gsl_matrix_get_col(res_column, matrix, n - 1);
      gsl_matrix_get_col(lead_column, matrix, lead_index_column);

      gsl_vector_div(res_column, lead_column);

      lead_index_row = 0;
      for (int i = 0; i < res_column->size; ++i)
      {
         double elem = gsl_vector_get(res_column, i);
         double min_elem = gsl_vector_get(res_column, lead_index_row);

         if (min_elem < 0 && elem > 0)
         {
            lead_index_row = i;
            continue;
         }

         if (elem < min_elem && elem > 0)
         {
            lead_index_row = i;
         }
      }

      lead_elem = gsl_matrix_get(matrix, lead_index_row, lead_index_column);
      if (lead_elem <= 0)
      {
         fprintf(stderr, PRINT_FORMAT_ERROR "Lead element should not be less than zero!\n" PRINT_CLEAR_FORMAT);
         goto done;
      }
      lead_row = gsl_matrix_row(matrix, lead_index_row);
      gsl_vector_scale(&lead_row.vector, 1.0 / lead_elem);

      for (int i = 0; i < m; ++i)
      {
         if (i == lead_index_row)
         {
            continue;
         }

         current_row = gsl_matrix_row(matrix, i);
         gsl_matrix_get_row(lead_row_copy, matrix, lead_index_row);

         gsl_vector_scale(lead_row_copy, gsl_vector_get(&current_row.vector,
                                                        lead_index_column));
         gsl_vector_sub(&current_row.vector, lead_row_copy);
      }

      dump_matrix_h_row(matrix, lead_index_row);
   }

   dump_matrix_result(matrix);

done:

   gsl_matrix_free(matrix);
   gsl_vector_free(lead_column);
   gsl_vector_free(res_column);
   gsl_vector_free(lead_row_copy);

   return GSL_SUCCESS;
}
