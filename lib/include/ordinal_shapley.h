#ifndef ORDINAL_SHAPLEY_H
#define ORDINAL_SHAPLEY_H

#include <vector>

#include "shapley.h"

namespace shapley
{

    class OrdinalPlayer : public Player
    {
    public:
        explicit OrdinalPlayer(int position_in, const std::vector<double> &v);
        int position() const;

    protected:
        int m_position;
        const std::vector<double> &v_ref;
    };

    class OrdinalCharacteristicFunction : public CharacteristicFunction<OrdinalPlayer>
    {
    public:
        explicit OrdinalCharacteristicFunction(const std::vector<double> &v);
        double value(const Coalition<OrdinalPlayer> &coalition) const override;

    protected:
        const std::vector<double> &v_ref;
    };

} // namespace shapley

#endif // ORDINAL_SHAPLEY_H
