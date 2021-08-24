// This file requires at least C99 to compile

/**
 * @file   facets.c
 * @author Jean-Cédric Chappelier <jean-cedric.chappelier@epfl.ch>
 * @author Merlin Nimier-David <merlin.nimier-david@epfl.ch>
 *
 * @copyright EPFL 2021
**/
/**
 * @section DESCRIPTION
 *
 * Template du homework du cours CS-207, année 2021.
**/

#include <limits.h>  // For INT_MAX
#include <math.h>
#include <stdio.h>
#include <stdlib.h>  // For EXIT_SUCCESS/FAILURE

// ----------------------------------------------
//   ___             _            _
//  / __|___ _ _  __| |_ __ _ _ _| |_ ___
// | (__/ _ \ ' \(_-<  _/ _` | ' \  _(_-<
//  \___\___/_||_/__/\__\__,_|_||_\__/__/

#define NB_VERTICES_MAX 128
#define NB_FACES_MAX    128
#define NB_FACETS_PER_VERTEX_MAX 20
#define NO_FACET ((Facet*) NULL)
#define NB_VERTICES_PER_FACE 3
#define NULL_VECTOR ((Vector) {0.0, 0.0, 0.0})

// ----------------------------------------------
//  _____
// |_   _|  _ _ __  ___ ___
//   | || || | '_ \/ -_|_-<
//   |_| \_, | .__/\___/__/
//       |__/|_|

typedef struct {
    double x;
    double y;
    double z;
} Vector;

typedef struct Facet Facet;

typedef struct {
    Vector sommet;
    Facet *facets[NB_FACETS_PER_VERTEX_MAX];
    size_t nb_facets;
    size_t nb_current_facets;
} Vertex;

struct Facet {
    Vertex *sommets[NB_VERTICES_PER_FACE];
};


Vector vertices_to_vector(Vertex a, Vertex b);

Vector vector_product(Vector a, Vector b);

Vector vector_norm(Vector a);

double norm(Vector a);

void vider_tampon(void);


// ======================================================================
/**
 * Attempts to register `p` as the `vertex_i`th vertex of face `f`,
 * and to add a reference to `f` in `p`.
 *
 * Returns 1 if successful, 0 otherwise.
 */
int register_vertex_for_face(Vertex *p, Facet *f, size_t vertex_i) {
    f->sommets[vertex_i] = p;
    p->facets[p->nb_current_facets] = f;

    if (f->sommets[vertex_i] != NULL && p->facets[p->nb_facets] != NULL) {
        return 1;
    }
    return 0;
}


// ======================================================================
/**
 * Lets the user input their desired number of 3D vertices (< nb_max)
 * by first giving the number of vertices and then giving the 3D coordinates
 * of each vertex.
 *
 * The vertices are stored in `tab`.
 * Returns the number of vertices succesfully received.
 */
size_t query_vertices(Vertex tab[], size_t nb_max) {
    size_t sommet = 0;

    do {
        fprintf(stdout, "Combien de sommets ?");
        fflush(stdout);
        int isConform = scanf("%zu", &sommet);
        if (isConform != 1) {
            putchar('\n');
            fprintf(stderr, "Erreur: entree un entier positif plus petit que %zu \n", nb_max);
            vider_tampon();
        }
    } while (!feof(stdin) && !ferror(stdin) && ((sommet < 1) || (sommet > nb_max)));

    putchar('\n');

    for (size_t i = 0; i < sommet; ++i) {
        fprintf(stdout, "    Vertex %zu ? ", i);
        int isXConform = scanf("%lf", &tab[i].sommet.x);
        if (isXConform != 1) {
            putchar('\n');
            fprintf(stderr, "Erreur: entree un normbre \n");
            vider_tampon();
            return i;
        }
        putchar(' ');
        int isYConform = scanf("%lf", &tab[i].sommet.y);
        if (isYConform != 1) {
            putchar('\n');
            fprintf(stderr, "Erreur: entree un normbre \n");
            vider_tampon();
            return i;
        }
        putchar(' ');
        int isZConform = scanf("%lf", &tab[i].sommet.z);
        if (isZConform != 1) {
            putchar('\n');
            fprintf(stderr, "Erreur: entree un normbre \n");
            vider_tampon();
            return i;
        }
        putchar('\n');
    }
    return sommet;
}

// ======================================================================
/**
 * Lets the user input their desired number of triangular
 * facets (< nb_max) by first giving the number of facets,
 * and then giving the index (0-based) of its 3 vertices.
 *
 * The facets are stored in `faces`, and `pts` is updated
 * with pointers to the facets referring to it.
 * Returns the number of facets successfully received.
 */
size_t query_facets(Facet faces[], size_t nb_max,
                    Vertex pts[], size_t nb_vertices) {
    size_t facet = 0;

    for (int i = 0; i < NB_FACES_MAX; ++i) {
        for (int j = 0; j < NB_VERTICES_PER_FACE; ++j) {
            faces[i].sommets[j] = NULL;
        }
    }
    for (int i = 0; i < nb_vertices; ++i) {
        for (int j = 0; j < NB_FACETS_PER_VERTEX_MAX; ++j) {
            pts[i].facets[j] = NULL;
        }
    }

    do {
        fprintf(stdout, "Combien de facettes ? ");
        fflush(stdout);
        int isConform = scanf("%zu", &facet);
        if (isConform != 1) {
            putchar('\n');
            fprintf(stderr, "Erreur: entree un entier positif plus petit que %zu \n", nb_max);
            vider_tampon();
        }
    } while (!feof(stdin) && !ferror(stdin) && ((facet < 1) || (facet > nb_max)));

    putchar('\n');

    int point = 0;

    for (size_t i = 0; i < facet; ++i) {
        fprintf(stdout, "    Face %zu ? ", i);
        for (size_t j = 0; j < NB_VERTICES_PER_FACE; ++j) {
            int isConform = scanf("%d", &point);
            if (point < 0 || point >= nb_vertices || isConform != 1) {
                putchar('\n');
                fprintf(stderr, "Erreur: les indices des sommets doivent etre entre 0 et %zu \n", nb_vertices);
                vider_tampon();
                return i;
            } else {
                pts[point].nb_current_facets = i;
                pts[point].nb_facets +=
                        register_vertex_for_face(&pts[point], &faces[i], j);
            }
        }
        putchar('\n');
    }
    return facet;
}

/**
 * Fonction qui vide le tampon
 */
void vider_tampon(void) {
    while (!feof(stdin) && !ferror(stdin) && getc(stdin) != '\n');
}

// ======================================================================
#define PREC 1e-10
#define round(x) (fabs(x) < PREC ? 0.0 : x)

void vector_println(const Vector *pv) {
    printf("(%.4g, %.4g, %.4g)\n", round(pv->x), round(pv->y), round(pv->z));
}

// ======================================================================
Vector facet_normal(size_t index, const Facet t[], size_t nb, double *out_area) {
    Vector rez = NULL_VECTOR;

    if (index >= nb) {
        fprintf(stderr,
                "Erreur: la valeur de l'index donné [%zu] est superieur au nombre d'éléments dans t [%zu]\n",
                index, nb);
        return rez;
    }

    Vector a = vertices_to_vector(*t[index].sommets[1],
                                  *t[index].sommets[0]);

    Vector b = vertices_to_vector(*t[index].sommets[2],
                                  *t[index].sommets[1]);

    Vector c = vertices_to_vector(*t[index].sommets[1],
                                  *t[index].sommets[2]);

    rez = vector_norm(vector_product(a, c));

    if (out_area != NULL) {
        *out_area = norm(vector_product(a, b)) / 2;
    }
    return rez;
}

/**
 * Calcule un vecteur pointant du sommet a au sommet b
 *
 */
Vector vertices_to_vector(const Vertex a, const Vertex b) {
    Vector ret = {a.sommet.x - b.sommet.x,
                  a.sommet.y - b.sommet.y,
                  a.sommet.z - b.sommet.z};
    return ret;
}

/**
 * Calcule le produit vectoriel entre a et b
 *
 */
Vector vector_product(const Vector a, const Vector b) {
    Vector ret = {a.y * b.z - a.z * b.y,
                  a.z * b.x - a.x * b.z,
                  a.x * b.y - a.y * b.x};
    return ret;
}

/**
 * Calcule le vecteur normale de a
 *
 */
Vector vector_norm(Vector a) {
    double norme = norm(a);
    if (norm(a) != 0) {
        a.x /= norme;
        a.y /= norme;
        a.z /= norme;
    }
    return a;
}

/**
 * Calcule la norme de a
 *
 */
double norm(const Vector a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

// ======================================================================
Vector mesh_barycenter(const Vertex vertices[], size_t nb_vertices) {
    Vector rez = NULL_VECTOR;
    if (nb_vertices > 0) {
        for (int i = 0; i < nb_vertices; ++i) {
            rez.x += vertices[i].sommet.x;
            rez.y += vertices[i].sommet.y;
            rez.z += vertices[i].sommet.z;
        }
        rez.x /= nb_vertices;
        rez.y /= nb_vertices;
        rez.z /= nb_vertices;
    }
    return rez;
}

// ======================================================================
void print_mesh_stats(const Vertex p[], size_t nb_vertices,
                      const Facet f[], size_t nb_facets) {
    const Vector bary = mesh_barycenter(p, nb_vertices);
    putchar('\n');
    printf("La forme contient %zu sommets et %zu facettes.\n", nb_vertices, nb_facets);
    printf("Le barycentre est : ");
    vector_println(&bary);
    putchar('\n');
}

// ======================================================================
void print_vertex_info(size_t index, const Vertex p[], size_t nb) {
    if (index >= nb) {
        fprintf(stderr,
                "Erreur: la valeur de l'index donné [%zu] est superieur au nombre d'éléments dans t [%zu]\n",
                index, nb);
        return;
    }
    printf("Le sommet %03zu appartient à %zu facettes : ", index, p[index].nb_facets);

    for (size_t i = 0; i < NB_FACETS_PER_VERTEX_MAX; ++i) {
        if (p[index].facets[i] != NO_FACET) {
            fprintf(stdout, "%03zu ", i);
        }
    }
    putchar('\n');
}

// ======================================================================
void print_face_info(size_t index, const Facet f[], size_t nb) {
    if (index >= nb) {
        fprintf(stderr, "Erreur : impossible de calculer la normale à la face %zu\n",
                index);
        fprintf(stderr, "Erreur : les indices des facettes sont entre 0 et %zu\n", nb - 1);
        return;
    }

    double area = 0.0;
    const Vector n = facet_normal(index, f, nb, &area);
    printf("Face %03zu : aire = %.4g, normale = ", index, area);
    vector_println(&n);
}

// ----------------------------------------------
//  __  __      _
// |  \/  |__ _(_)_ _
// | |\/| / _` | | ' \
// |_|  |_\__,_|_|_||_|

#define MAX_PRINT 10

int main(void) {
    Vertex vertices[NB_VERTICES_MAX];
    const size_t nb_vertices = query_vertices(vertices, NB_VERTICES_MAX);

    Facet facets[NB_FACES_MAX];
    const size_t nb_facets = query_facets(facets, NB_FACES_MAX,
                                          vertices, nb_vertices);

    print_mesh_stats(vertices, nb_vertices, facets, nb_facets);

    for (size_t i = 0; i < (nb_vertices < MAX_PRINT ? nb_vertices : MAX_PRINT); ++i) {
        print_vertex_info(i, vertices, nb_vertices);
    }

    for (size_t i = 0; i < nb_facets; ++i) {
        print_face_info(i, facets, nb_facets);
    }

    return EXIT_SUCCESS;
}
