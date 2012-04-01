/**
 * @file
 * @author Alexander Sherikov
 * @date 30.03.2012 14:35:57 MSD
 */


#ifndef CONSTRAINT_H
#define CONSTRAINT_H

/****************************************
 * INCLUDES 
 ****************************************/



/****************************************
 * TYPEDEFS 
 ****************************************/
/// @addtogroup gAS
/// @{

class active_constraint
{
    public:
        active_constraint (int ind_, int sign_)
        {
            ind = ind_;
            sign = sign_;
        }

        int ind;
        int sign;
};


/** \brief Defines constraints associated with states of the system. */
class constraint
{
    public:
        /**
         * @brief Set parameters of the bound
         *
         * @param[in] ind_ index of the first element of a state in the vector of states
         * @param[in] coef_x_ coefficient for x coordinate
         * @param[in] coef_y_ coefficient for y coordinate
         * @param[in] lb_ lower bound
         * @param[in] ub_ upper bound
         * @param[in] active activity of the bound
         */
        void set(
                const int ind_, 
                const double coef_x_, 
                const double coef_y_, 
                const double lb_, 
                const double ub_, 
                const bool active)
        {
            ind = ind_;
            coef_x = coef_x_;
            coef_y = coef_y_;
            lb = lb_;
            ub = ub_;
            isActive = active;
        }


        /** Variable number (on which to impose the bounds). */
        int ind;

        //@{
        /// Coefficients
        double coef_x;
        double coef_y;
        //@}


        /** Lower bound. */
        double lb;

        /** Upper bound. */
        double ub;

        /** 
         * If isActive then one of the bounds is in the 
         * working set.
         */
        bool isActive;
};


///@}

#endif /*CONSTRAINT_H*/

