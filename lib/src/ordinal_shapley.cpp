#include "ordinal_shapley.h"

using namespace std;

namespace Shapley
{

    OrdinalPlayer::OrdinalPlayer(int position_in, const vector<double> &v)
        : position(position_in), v_ref(v) {}

    int OrdinalPlayer::getPosition() const
    {
        return position;
    }

    OrdinalCharacteristicFunction::OrdinalCharacteristicFunction(const vector<double> &v)
        : v_ref(v) {}

    double OrdinalCharacteristicFunction::getValue(const Coalition<OrdinalPlayer> &coalition) const
    {        
        int coalition_position = 0;
        for (const OrdinalPlayer *member : coalition.getMembers())
        {
            const int player_position = member->getPosition();
            coalition_position = coalition_position | (1 << player_position);
        }
        const double contribution = v_ref[coalition_position - 1];
        return contribution;
    }
}