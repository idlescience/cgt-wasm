#ifndef ORDINAL_SHAPLEY_H
#define ORDINAL_SHAPLEY_H

#include <vector>

#include "shapley.h"

namespace Shapley
{

    class OrdinalPlayer : public Player
    {
    public:
        explicit OrdinalPlayer(int position_in, const std::vector<double> &v);
        double getContribution() const override;
        int getPosition() const;

    protected:
        int position;
        const std::vector<double> &v_ref;
    };

    class OrdinalCharacteristicFunction : public CharacteristicFunction<OrdinalPlayer>
    {
    public:
        explicit OrdinalCharacteristicFunction(const std::vector<double> &v);
        double getValue(const Coalition<OrdinalPlayer> &coalition) const override;

    protected:
        const std::vector<double> &v_ref;
    };

}

#endif // ORDINAL_SHAPLEY_H
