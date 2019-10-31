#ifndef MACCSKEYS_FN2_H
#define MACCSKEYS_FN2_H

#define NUMBER_MACCS_KEYS 192

class MK_Molecular_Properties;
class MACCSKeys
{
  private:

    int _aromatic_nitrogens_do_not_have_hydrogens;

    int _nbits;  

//  private functions

    void _default_values();

    int  _key0(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key1(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key2(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key3(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key4(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key5(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key6(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key7(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key8(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key9(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key10(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key11(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key12(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key13(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key14(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key15(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key16(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key17(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key18(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key19(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key20(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key21(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key22(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key23(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key24(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key25(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key26(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key27(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key28(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key29(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key30(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key31(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key32(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key33(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key34(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key35(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key36(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key37(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key38(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key39(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key40(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key41(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key42(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key43(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key44(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key45(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key46(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key47(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key48(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key49(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key50(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key51(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key52(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key53(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key54(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key55(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key56(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key57(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key58(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key59(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key60(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key61(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key62(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key63(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key64(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key65(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key66(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key67(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key68(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key69(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key70(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key71(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key72(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key73(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key74(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key75(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key76(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key77(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key78(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key79(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key80(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key81(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key82(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key83(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key84(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key85(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key86(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key87(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key88(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key89(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key90(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key91(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key92(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key93(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key94(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key95(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key96(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key97(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key98(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key99(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key100(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key101(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key102(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key103(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key104(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key105(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key106(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key107(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key108(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key109(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key110(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key111(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key112(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key113(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key114(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key115(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key116(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key117(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key118(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key119(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key120(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key121(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key122(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key123(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key124(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key125(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key126(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key127(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key128(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key129(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key130(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key131(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key132(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key133(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key134(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key135(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key136(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key137(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key138(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key139(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key140(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key141(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key142(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key143(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key144(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key145(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key146(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key147(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key148(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key149(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key150(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key151(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key152(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key153(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key154(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key155(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key156(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key157(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key158(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key159(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key160(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key161(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key162(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key163(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key164(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key165(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key166(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key167(const Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key168(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key169(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key170(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key171(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key172(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key173(Molecule & m, const MK_Molecular_Properties & mpr) const;
    void _key174(Molecule & m, const MK_Molecular_Properties & mpr, int *) const;
    int  _key175(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key176(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key177(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key178(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key179(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key180(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key181(Molecule & m, const MK_Molecular_Properties & mpr) const;
    void _key182(Molecule & m, const MK_Molecular_Properties & mpr, int *) const;
    int  _key184(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key185(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key186(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key187(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key188(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key189(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key190(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key191(Molecule & m, const MK_Molecular_Properties & mpr) const;
    int  _key192(Molecule & m, const MK_Molecular_Properties & mpr) const;

    int _is_para_substituted_within_scaffold (const Molecule & m,
                                             const Set_of_Atoms & r,
                                             const MK_Molecular_Properties & mpr) const;

  public:
    MACCSKeys ();

    int nbits () const { return _nbits;}

    int set_nbits (int);

    void set_aromatic_nitrogens_do_not_have_hydrogens(int s) { _aromatic_nitrogens_do_not_have_hydrogens = s;}

//  someone has already generated a fingerprint and wants the level 2 variant

    int set_level_2_fingerprint (int *) const;

    int operator () (Molecule &, int *) const;
};

#endif
