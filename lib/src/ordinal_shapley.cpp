#include "ordinal_shapley.h"

using namespace std;

namespace shapley
{

    OrdinalPlayer::OrdinalPlayer(int position_in)
        : m_position(position_in) {}

    int OrdinalPlayer::position() const
    {
        return m_position;
    }

    OrdinalCharacteristicFunction::OrdinalCharacteristicFunction(const vector<double> &v)
        : v_ref(v) {}

    double OrdinalCharacteristicFunction::value(const Coalition<OrdinalPlayer> &coalition) const
    {        
        int coalition_position = 0;
        for (const OrdinalPlayer *member : coalition.getMembers())
        {
            const int player_position = member->position();
            coalition_position = coalition_position | (1 << player_position);
        }
        const double contribution = v_ref[coalition_position - 1];
        return contribution;
    }
}