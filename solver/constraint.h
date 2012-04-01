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

/** \brief Defines constraints associated with states of the system. */
class constraint
{
    public:
        /**
         * @brief Set parameters of the bound
         *
         * @param[in] cind_ the number of constraint
         * @param[in] coef_x_ coefficient for x coordinate
         * @param[in] coef_y_ coefficient for y coordinate
         * @param[in] lb_ lower bound
         * @param[in] ub_ upper bound
         * @param[in] active activity of the bound
         */
        void set(
                const int cind_, 
                const double coef_x_, 
                const double coef_y_, 
                const double lb_, 
                const double ub_, 
                const bool active)
        {
            cind = cind_;
            ind = cind/2*SMPC_NUM_STATE_VAR;
            coef_x = coef_x_;
            coef_y = coef_y_;
            lb = lb_;
            ub = ub_;
            isActive = active;
        }


        /// The sequential number of the constraint
        int cind;

        /// Index of the first element of a state in the vector of states
        int ind;


        //@{
        /// Coefficients
        double coef_x;
        double coef_y;
        //@}


        /// Lower bound.
        double lb;

        /// Upper bound.
        double ub;


        /** 
         * Since we do not distinguish lower/upper bounds of active constraints (<= and => 
         * inequlities are treated in the same way), we have to adjust signs of Lagrange 
         * multipliers before downdate. See also '@ref pBounds'.
         */
        int sign;

        /// If isActive then one of the bounds is in the working set.
        bool isActive;
};


///@}

#endif /*CONSTRAINT_H*/

