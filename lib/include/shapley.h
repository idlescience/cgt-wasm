#ifndef SHAPLEY_LIBRARY_H
#define SHAPLEY_LIBRARY_H

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

#include "group.h"

using namespace Group;

/**
 * Implementation of the computation of the game theoretic Shapley value.
 */
namespace Shapley
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
            bool *contain = new bool[n]{0};

            for (int i = 0; i < n; i++)
            {
                contain[i] = 1;
                do
                {
                    Coalition<PlayerType> input;
                    for (int j = 0; j < n; j++)
                    {
                        if (contain[j])
                        {
                            input.add(members.at(j));
                        }
                    }
                    ans.push_back(input);
                } while (std::prev_permutation(contain, contain + n));
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
    static std::map<const PlayerType *, double> compute(const std::vector<const PlayerType *> &players,
                                                        const CharacteristicFunction<PlayerType> &char_func)
    {
        const int n = players.size();
        std::map<const PlayerType *, double> shapley_values;
        for (size_t i = 0; i < n; i++)
        {
            shapley_values[players.at(i)] = 0.0;
        }

        Coalition<PlayerType> grand_coalition(players);
        std::vector<Coalition<PlayerType>> B = grand_coalition.power_set();

        for (Coalition<PlayerType> &coalition : B)
        {
            for (size_t i = 0; i < coalition.size(); i++)
            {
                const PlayerType *member = coalition.player_at(i);
                if (coalition.size() > 1)
                {
                    shapley_values[member] += (1 / binomial_coef(n, coalition.size())) *
                                              (char_func.value(coalition) - char_func.value(coalition.except(i)));
                }
                else
                {
                    shapley_values[member] += char_func.value(coalition);
                }
            }
        }

        for (size_t i = 0; i < n; i++)
        {
            shapley_values[players.at(i)] /= n;
        }

        return shapley_values;
    }
} // namespace Shapley

#endif
