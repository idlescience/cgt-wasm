#ifndef SHAPLEY_LIBRARY_H
#define SHAPLEY_LIBRARY_H

#include <algorithm>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

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

        /**
         * @return This player's contribution to the characteristic function.
         */
        virtual double getContribution() const = 0;
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
        Coalition getUpUntil(size_t index)
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
        Coalition getExcept(size_t index)
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
        const PlayerType *getPlayerAt(size_t index)
        {
            return members.at(index);
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

        virtual double getValue(const Coalition<PlayerType> &coalition) const = 0;
    };

    /** IMPLEMENTATION */

    static unsigned nChoosek(unsigned n, unsigned k)
    {
        if (k > n)
        {
            return 0;
        }
        if (k * 2 > n)
        {
            k = n - k;
        }
        if (k == 0)
        {
            return 1;
        }

        int result = n;
        for (int i = 2; i <= k; ++i)
        {
            result *= (n - i + 1);
            result /= i;
        }
        return result;
    }

    /**
     * Computes each player's Shapley value.
     * @tparam PlayerType A class that derives from Shapley::Player.
     * @param players
     * @param charFunc A characteristic function for 'PlayerType'.
     * @return A [Player, Shapley Value]-map.
     */
    template <class PlayerType>
    static std::map<const PlayerType *, double> compute(const std::vector<const PlayerType *> &players,
                                                        const CharacteristicFunction<PlayerType> &charFunc)
    {
        const int n = players.size();
        std::map<const PlayerType *, double> shapley_values;
        for (size_t i = 0; i < n; i++)
        {
            shapley_values[players.at(i)] = 0.0;
        }

        Coalition<PlayerType> grand_coalition(players);
        std::vector<Coalition<PlayerType>> power_set = powerSet(grand_coalition);

        for (Coalition<PlayerType> &coalition : power_set)
        {
            for (size_t i = 0; i < coalition.size(); i++)
            {
                const PlayerType *member = coalition.getPlayerAt(i);
                if (coalition.size() > 1)
                {
                    shapley_values[member] +=
                        (1 / nChoosek(n, coalition.size())) *
                        (charFunc.getValue(coalition) - charFunc.getValue(coalition.getExcept(i)));
                }
                else
                {
                    shapley_values[member] += charFunc.getValue(coalition);
                }
            }
        }

        for (size_t i = 0; i < n; i++)
        {
            shapley_values[players.at(i)] /= n;
        }

        return shapley_values;
    }

    /**
     * Retrieves powerset of a coalition.
     * @tparam PlayerType A class that derives from Shapley::Player.
     * @param coalition
     * @return A PlayerType std::vector.
     */
    template <class PlayerType> static std::vector<Coalition<PlayerType>> powerSet(Coalition<PlayerType> &coalition)
    {
        int n = coalition.size();
        std::vector<Coalition<PlayerType>> ans = {};

        for (int i = 0; i < n; i++)
        {
            const PlayerType *element = coalition.getPlayerAt(i);
            int len = ans.size();

            for (int j = 0; j < len; j++)
            {
                Coalition<PlayerType> temp = ans[j];
                temp.add(element);
                ans.push_back(temp);
            }
        }

        return ans;
    }
} // namespace Shapley

#endif
