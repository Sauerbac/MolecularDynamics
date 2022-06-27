#ifndef NEIGHBORS_H
#define NEIGHBORS_H

#include "Atoms.h"

class NeighborList {
  public:
    NeighborList(double interaction_range);

    /*
     * Update neighbor list from the particle positons stores in the `atoms`
     * argument
     */
    const std::tuple<const Eigen::ArrayXi &, const Eigen::ArrayXi &>
    update(const Atoms &atoms);

    /*
     * Return internal seed and neighbor arrays
     */
    const std::tuple<const Eigen::ArrayXi &, const Eigen::ArrayXi &>
    neighbors() const {
        if (seed_.size() > 0) {
            return {seed_, neighbors_};
        } else {
            throw std::runtime_error(
                "Neighbor list not yet computed. Use `update` to compute it.");
        }
    }

    /*
     * Return the total number of neighbors found by the last call to `update`
     */
    int nb_neighbors() const {
        // Note that the length of the `neighbor_` array can be longer than the
        // total number of neighbors!
        return seed_(seed_.size() - 1);
    }

    double get_interaction_range() {
        return interaction_range_;
    }

    /*
     * Return the number of neighbors of atom `i` found by the last call to
     * `update`
     */
    int nb_neighbors(int i) const {
        assert(i >= 0);
        assert(i < seed_.size());
        return seed_(i + 1) - seed_(i);
    }

    class iterator
        : public std::iterator<std::input_iterator_tag,
                               std::tuple<int, int>,   // value_type
                               std::tuple<int, int>,   // difference_type
                               std::tuple<int, int> *, // pointer
                               std::tuple<int, int>> { // reference
      public:
        explicit iterator(const Eigen::ArrayXi &seed,
                          const Eigen::ArrayXi &neighbors, int i, int n)
            : seed_{seed},
              neighbors_{neighbors},
              i_{i},
              n_{n} {
        }

        iterator &operator++() {
            assert(n_ < seed_(seed_.size() - 1));

            n_++;
            while (n_ == seed_(i_ + 1) && i_ < seed_.size() - 2) {
                i_++;
            }
            return *this;
        }

        bool operator==(iterator other) const {
            return n_ == other.n_;
        }

        bool operator!=(iterator other) const {
            return n_ != other.n_;
        }

        reference operator*() const {
            return {i_, neighbors_(n_)};
        }

      protected:
        const Eigen::ArrayXi &seed_;
        const Eigen::ArrayXi &neighbors_;
        int i_, n_;
    };

    /*
     * Return iterator that represents the beginning of the neighbor list
     */
    iterator begin() const {
        if (seed_.size() > 0) {
            int i0{0}, n0{seed_(0)};
            assert(n0 == 0); // seed_(0) should be zero
            // If the first atom has no neighbors, we need to fast-forward the
            // counter to the first atom that does have a neighbor.
            while (n0 == seed_(i0 + 1) && n0 < nb_neighbors()) {
                i0++;
            }
            return iterator(seed_, neighbors_, i0, n0);
        } else {
            throw std::runtime_error(
                "Neighbor list not yet computed. Use `update` to compute it.");
        }
    }

    /*
     * Return iterator that represents the end of the neighbor list
     */
    iterator end() const {
        if (seed_.size() > 0) {
            return iterator(seed_, neighbors_, seed_.size() - 2,
                            nb_neighbors());
        } else {
            throw std::runtime_error(
                "Neighbor list not yet computed. Use `update` to compute it.");
        }
    }

  protected:
    template <typename T>
    static decltype(auto)
    coordinate_to_index(const T &x, const T &y, const T &z,
                        const Eigen::Array3i &nb_grid_pts) {
        return x + nb_grid_pts(0) * (y + nb_grid_pts(1) * z);
    }

    int coordinate_to_index(const Eigen::Array3i &c,
                            const Eigen::Array3i &nb_grid_pts) {
        return coordinate_to_index(c(0), c(1), c(2), nb_grid_pts);
    }

    double interaction_range_;

    Eigen::ArrayXi seed_;
    Eigen::ArrayXi neighbors_;
};

#endif // NEIGHBORS_H