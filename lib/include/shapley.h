#ifndef SHAPLEY_LIBRARY_H
#define SHAPLEY_LIBRARY_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

#include "group.h"

using namespace group;

/**
 * Implementation of the computation of the game theoretic Shapley value.
 */
namespace shapley
{

    /** REQUIRED CLASSES */

    /**
     * A Player contributes to a coalition's worth in some way.
     */
    class Player
    {
    public:
        virtual ~Player()
        {
        }
    };

    /**
     * A Coalition contains any number of Players.
     */
    template <class PlayerType> class Coalition
    {
    public:
        /**
         * Constructs an empty coalition.
         */
        Coalition<PlayerType>()
        {
        }

        /**
         * @param members This coalition's members.
         * The Coalition will *not* delete the pointers upon destruction.
         */
        explicit Coalition<PlayerType>(const std::vector<const PlayerType *> &members)
        {
            // Call copy constructor.
            this->members = std::vector<const PlayerType *>(members);
        }

        virtual ~Coalition<PlayerType>()
        {
        }

        /**
         * Adds 'member' to this coalition.
         * @param member
         * @throws invalid_argument If the coalition already contains 'member'.
         */
        void add(const PlayerType *member)
        {
            if (contains(member))
            {
                throw std::invalid_argument("Coalition::add called, but 'member' is already contained.");
            }
            this->members.push_back(member);
        }

        void remove(const PlayerType *member)
        {
            members.erase(std::remove(members.begin(), members.end(), member), members.end());
        }

        bool contains(const PlayerType *member) const
        {
            return std::find(members.begin(), members.end(), member) != members.end();
        }

        size_t size() const
        {
            return members.size();
        }

        const std::vector<const PlayerType *> &getMembers() const
        {
            return members;
        }

        void clear()
        {
            members.clear();
        }

        /**
         * @param index
         * @return A Coalition with all members up to but not including 'index'.
         */
        Coalition until(size_t index)
        {
            Coalition copy;
            for (size_t i = 0; i < index; i++)
            {
                copy.add(members.at(i));
            }
            return copy;
        }

        /**
         * @param index
         * @return A Coalition with all members not including 'index'.
         */
        Coalition except(size_t index)
        {
            Coalition copy;
            for (size_t i = 0; i < members.size(); i++)
            {
                if (i != index)
                {
                    copy.add(members.at(i));
                }
            }
            return copy;
        }

        /**
         * @param index
         * @return the Player at 'index'.
         */
        const PlayerType *player_at(size_t index)
        {
            return members.at(index);
        }

        /**
         * Retrieves powerset of the coalition.
         * @return A PlayerType std::vector.
         */
        const std::vector<Coalition<PlayerType>> power_set()
        {
            const int n = members.size();
            std::vector<Coalition<PlayerType>> ans = {};
            bool *contains = new bool[n]{0};

            for (int i = 0; i < n; i++)
            {
                contains[i] = 1;
                do
                {
                    Coalition<PlayerType> input;
                    for (int j = 0; j < n; j++)
                    {
                        if (contains[j])
                        {
                            input.add(members.at(j));
                        }
                    }
                    ans.push_back(input);
                } while (std::prev_permutation(contains, contains + n));
            }
            return ans;
        }

    protected:
        std::vector<const PlayerType *> members;
    };

    /**
     * The Characteristic Function determines a Coalition's worth.
     */
    template <class PlayerType> class CharacteristicFunction
    {
    public:
        virtual ~CharacteristicFunction()
        {
        }

        virtual double value(const Coalition<PlayerType> &coalition) const = 0;
    };

    /**
     * Computes each player's Shapley value.
     * @tparam PlayerType A class that derives from Shapley::Player.
     * @param players
     * @param char_func A characteristic function for 'PlayerType'.
     * @return A [Player, Shapley Value]-map.
     */
    template <class PlayerType>
    static std::map<unsigned short int, double> compute(const std::vector<const PlayerType *> &players,
                                                        const CharacteristicFunction<PlayerType> &char_func)
    {
        Coalition<PlayerType> grand_coalition(players);
        const int n = grand_coalition.size();

        std::map<const PlayerType *, double> shapley_values_refs;
        for (unsigned short int i = 0; i < n; i++)
        {
            shapley_values_refs[grand_coalition.player_at(i)] = 0.0;
        }

        std::vector<Coalition<PlayerType>> B = grand_coalition.power_set();

        for (Coalition<PlayerType> &coalition : B)
        {
            const unsigned short int s = coalition.size();
            for (size_t i = 0; i < s; i++)
            {
                const PlayerType *member = coalition.player_at(i);
                const double full_value = char_func.value(coalition);
                const double binomial = binomial_coef(n - 1, s - 1);
                const double inv_binomial = 1 / binomial;
                double addend = 0;

                if (s > 1)
                {
                    const double skipped_value = char_func.value(coalition.except(i));
                    const double diff = full_value - skipped_value;
                    addend = inv_binomial * diff;

                    if (std::isinf(addend))
                    {
                        throw std::runtime_error("Shapley value is infinite.");
                    }
                }
                else
                {
                    addend = inv_binomial * full_value;
                }

                shapley_values_refs[member] += addend;
            }
        }

        std::map<unsigned short int, double> shapley_values;
        for (unsigned short int i = 0; i < n; i++)
        {
            const double gross_value = shapley_values_refs[grand_coalition.player_at(i)];
            shapley_values[i] = gross_value / n;
        }

        return shapley_values;
    }
} // namespace shapley

#endif
