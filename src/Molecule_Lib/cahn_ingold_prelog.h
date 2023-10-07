  private:

    CahnIngoldPrelog CahnIngoldPrelogValue(const Chiral_Centre* c,
                                int top_front, int top_back, 
                                int left_down, int right_down);
    CahnIngoldPrelog CahnIngoldPrelogValue(int north, int se, int sw) const;

    std::optional<uint32_t> ChiralCentreMemberToCipInt(int zatom) const;
